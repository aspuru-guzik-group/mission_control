import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings

from mc_utils import test_utils
from missions.models import Flow as FlowModel
from flow_engines.flow import Flow
from flow_runners.base_flow_runner import BaseFlowRunner
from flow_engines.flow_engine import FlowEngine
from flow_client.flow_client import MissionControlFlowClient


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('missions.urls')),
]

@override_settings(ROOT_URLCONF=__name__)
class FlowRunnerE2ETestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client_to_use_json(client=self.client)
        self.keyed_task_classes = self.generate_keyed_task_classes()
        self.flow_engine = self.generate_flow_engine(
            keyed_task_classes=self.keyed_task_classes)
        self.flow_client = MissionControlFlowClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )
        self.flow_runner = self.generate_flow_runner(
            flow_client=self.flow_client,
            flow_engine=self.flow_engine)
        self.populate_flows()

    def generate_keyed_task_classes(self):
        class BaseTestTask():
            def __init__(self, *args, **kwargs):
                self.__dict__.update(kwargs)

            def tick(self, *args, **kwargs):
                if not hasattr(self.__class__, 'ticks'):
                    self.__class__.ticks = 0
                self.__class__.ticks += 1

        task_classes = []
        class Task_0(BaseTestTask):
            def tick(self, *args, **kwargs):
                super().tick()
                next_task = Task_1(status='PENDING')
                self.flow.add_task(next_task, precursor=self)
                self.status = 'COMPLETED'
        task_classes.append(Task_0)

        class Task_1(BaseTestTask):
            allow_completion = False

            def tick(self, *args, **kwargs):
                super().tick()
                if self.__class__.allow_completion:
                    next_task = Task_2(status='PENDING')
                    self.flow.add_task(next_task, precursor=self)
                    self.status = 'COMPLETED'
        task_classes.append(Task_1)

        class Task_2(BaseTestTask):
            allow_completion = False

            def tick(self, *args, **kwargs):
                super().tick()
                if self.__class__.allow_completion:
                    self.status = 'COMPLETED'
        task_classes.append(Task_2)

        keyed_task_classes = {task_class.__name__: task_class
                              for task_class in task_classes}
        return keyed_task_classes

    def generate_flow_engine(self, keyed_task_classes=None):
        flow_engine = FlowEngine()
        for task_class in keyed_task_classes.values():
            flow_engine.register_task_class(task_class=task_class)
        return flow_engine

    def generate_flow_runner(self, flow_engine=None,
                                 flow_client=None):
        flow_runner = BaseFlowRunner(flow_engine=flow_engine,
                                             flow_client=flow_client)
        return flow_runner

    def populate_flows(self):
        flow = Flow()
        root_task = self.keyed_task_classes['Task_0'](flow=flow,
                                                      status='PENDING')
        flow.add_task(task=root_task, as_root=True)
        serialization = self.flow_engine.serialize_flow(flow)
        flow_model = FlowModel.objects.create(
            serialization=json.dumps(serialization))
        self.flow_uuid = flow_model.uuid

    def test_run_flows(self):
        states = []
        states.append(self.get_state(prev_state=None))

        def tick_and_capture_state(label=None):
            #print("tick_and_capture", label)
            self.flow_runner.tick()
            states.append(self.get_state(prev_state=states[-1]))

        tick_and_capture_state(label='a')
        self.assertEqual(
            states[-1]['task_classes']['Task_0']['ticks'], 1)

        tick_and_capture_state(label='b')
        self.assertEqual(
            states[-1]['task_classes']['Task_0']['ticks'], 1)
        self.assertEqual(
            states[-1]['task_classes']['Task_1']['ticks'], 1)

        tick_and_capture_state(label='c')
        self.assertEqual(
            states[-1]['task_classes']['Task_1']['ticks'], 2)
        self.assertFalse(
            'ticks' in states[-1]['task_classes']['Task_2'])

        self.keyed_task_classes['Task_1'].allow_completion = True
        tick_and_capture_state(label='d')
        self.assertEqual(
            states[-1]['task_classes']['Task_1']['ticks'], 3)

        tick_and_capture_state(label='e')
        self.assertEqual(
            states[-1]['task_classes']['Task_1']['ticks'], 3)
        self.assertEqual(
            states[-1]['task_classes']['Task_2']['ticks'], 1)

        self.keyed_task_classes['Task_2'].allow_completion = True
        tick_and_capture_state(label='f')
        self.assertEqual(
            states[-1]['task_classes']['Task_2']['ticks'], 2)
        self.assertEqual(
            states[-1]['flows'][self.flow_uuid]['status'], 'COMPLETED')

        tick_and_capture_state(label='g')
        self.assertEqual(
            states[-2]['flows'][self.flow_uuid]['modified'],
            states[-1]['flows'][self.flow_uuid]['modified'])

    def get_state(self, prev_state=None):
        state = {
            'flows': {
                flow_model.uuid: flow_model.__dict__
                for flow_model in FlowModel.objects.all()
            },
            'task_classes': {
                key: task_class.__dict__
                for key, task_class in self.keyed_task_classes.items()
            }
        }
        return state
