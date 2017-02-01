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
        test_utils.patch_request_client(request_client=self.client)
        self.keyed_task_engines = self.generate_keyed_task_engines()
        self.flow_engine = self.generate_flow_engine(
            keyed_task_engines=self.keyed_task_engines)
        self.flow_client = MissionControlFlowClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )
        self.flow_runner = self.generate_flow_runner(
            flow_client=self.flow_client,
            flow_engine=self.flow_engine)
        self.populate_flows()

    def generate_keyed_task_engines(self):
        class BaseTaskEngine():
            def __init__(self, *args, **kwargs):
                self.__dict__.update(kwargs)

            @classmethod
            def tick_task(cls, *args, **kwargs):
                if not hasattr(cls, 'ticks'):
                    cls.ticks = 0
                cls.ticks += 1

        task_engines = []
        class TaskEngine_0(BaseTaskEngine):
            @classmethod
            def tick_task(cls, *args, task=None, ctx=None, **kwargs):
                super().tick_task(*args, task=task, **kwargs)
                next_task = {'task_engine': TaskEngine_1.__name__,
                             'status': 'PENDING'}
                ctx['flow'].add_task(next_task, precursor_keys=[task['key']])
                task['status'] = 'COMPLETED'
        task_engines.append(TaskEngine_0)

        class TaskEngine_1(BaseTaskEngine):
            allow_completion = False

            @classmethod
            def tick_task(cls, *args, task=None, ctx=None, **kwargs):
                super().tick_task(*args, task=task, **kwargs)
                if cls.allow_completion:
                    next_task = {'task_engine': TaskEngine_2.__name__,
                                 'status': 'PENDING'}
                    ctx['flow'].add_task(next_task,
                                         precursor_keys=[task['key']])
                    task['status'] = 'COMPLETED'
        task_engines.append(TaskEngine_1)

        class TaskEngine_2(BaseTaskEngine):
            allow_completion = False

            @classmethod
            def tick_task(cls, *args, task=None, ctx=None, **kwargs):
                super().tick_task(*args, task=task, **kwargs)
                if cls.allow_completion: task['status'] = 'COMPLETED'
        task_engines.append(TaskEngine_2)

        keyed_task_engines = {task_engine.__name__: task_engine
                              for task_engine in task_engines}
        return keyed_task_engines

    def generate_flow_engine(self, keyed_task_engines=None):
        flow_engine = FlowEngine()
        for task_engine in keyed_task_engines.values():
            flow_engine.register_task_engine(task_engine=task_engine)
        return flow_engine

    def generate_flow_runner(self, flow_engine=None, flow_client=None):
        flow_runner = BaseFlowRunner(flow_engine=flow_engine,
                                     flow_client=flow_client)
        return flow_runner

    def populate_flows(self):
        flow = Flow()
        root_task = {'task_engine': 'TaskEngine_0', 'status': 'PENDING'}
        flow.add_task(task=root_task, as_root=True)
        serialization = self.flow_engine.serialize_flow(flow)
        flow_model = FlowModel.objects.create(
            serialization=json.dumps(serialization))
        self.flow_uuid = flow_model.uuid

    def test_run_flows(self):
        states = []
        states.append(self.get_state(prev_state=None))

        def tick_and_capture_state(label=None):
            self.flow_runner.tick()
            states.append(self.get_state(prev_state=states[-1]))

        tick_and_capture_state(label='a')
        self.assertEqual(states[-1]['task_engines']['TaskEngine_0']['ticks'], 1)

        tick_and_capture_state(label='b')
        self.assertEqual(states[-1]['task_engines']['TaskEngine_0']['ticks'], 1)
        self.assertEqual(states[-1]['task_engines']['TaskEngine_1']['ticks'], 1)

        tick_and_capture_state(label='c')
        self.assertEqual(states[-1]['task_engines']['TaskEngine_1']['ticks'], 2)
        self.assertFalse('ticks' in states[-1]['task_engines']['TaskEngine_2'])

        self.keyed_task_engines['TaskEngine_1'].allow_completion = True
        tick_and_capture_state(label='d')
        self.assertEqual(states[-1]['task_engines']['TaskEngine_1']['ticks'], 3)

        tick_and_capture_state(label='e')
        self.assertEqual(states[-1]['task_engines']['TaskEngine_1']['ticks'], 3)
        self.assertEqual(states[-1]['task_engines']['TaskEngine_2']['ticks'], 1)

        self.keyed_task_engines['TaskEngine_2'].allow_completion = True
        tick_and_capture_state(label='f')
        self.assertEqual(states[-1]['task_engines']['TaskEngine_2']['ticks'], 2)
        self.assertEqual(
            states[-1]['flows'][self.flow_uuid]['status'], 'COMPLETED')

        tick_and_capture_state(label='g')
        self.assertEqual(states[-2]['flows'][self.flow_uuid]['modified'],
                         states[-1]['flows'][self.flow_uuid]['modified'])

    def get_state(self, prev_state=None):
        state = {
            'flows': {
                flow_model.uuid: flow_model.__dict__
                for flow_model in FlowModel.objects.all()
            },
            'task_engines': {
                key: task_engine.__dict__
                for key, task_engine in self.keyed_task_engines.items()
            }
        }
        return state
