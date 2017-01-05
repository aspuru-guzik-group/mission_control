import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings

from missions.models import Flow as FlowModel
from flow_engines.flow import Flow, BaseNode
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
        self.configure_request_client()
        self.keyed_node_classes = self.generate_keyed_node_classes()
        self.flow_engine = self.generate_flow_engine(
            keyed_node_classes=self.keyed_node_classes)
        self.flow_client = MissionControlFlowClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )
        self.flow_runner = self.generate_flow_runner(
            flow_client=self.flow_client,
            flow_engine=self.flow_engine)
        self.populate_flows()

    def configure_request_client(self):
        orig_patch = self.client.patch
        def json_patch(path, data=None, **kwargs):
            return orig_patch(path, json.dumps(data),
                              content_type='application/json', **kwargs)
        self.client.patch = json_patch

    def generate_keyed_node_classes(self):
        class BaseTestNode(BaseNode):
            def tick(self, *args, **kwargs):
                print("tick", self.__class__.__name__)
                if not hasattr(self.__class__, 'ticks'):
                    self.__class__.ticks = 0
                self.__class__.ticks += 1

        node_classes = []
        class Node_0(BaseTestNode):
            def tick(self, *args, **kwargs):
                super().tick()
                next_node = Node_1(status='PENDING')
                self.flow.add_node(next_node, precursors=[self])
                self.status = 'COMPLETED'
        node_classes.append(Node_0)

        class Node_1(BaseTestNode):
            allow_completion = False

            def tick(self, *args, **kwargs):
                super().tick()
                if self.__class__.allow_completion:
                    next_node = Node_2(status='PENDING')
                    self.flow.add_node(next_node, precursors=[self])
                    self.status = 'COMPLETED'
        node_classes.append(Node_1)

        class Node_2(BaseTestNode):
            allow_completion = False

            def tick(self, *args, **kwargs):
                super().tick()
                if self.__class__.allow_completion:
                    self.status = 'COMPLETED'
        node_classes.append(Node_2)

        keyed_node_classes = {node_class.__name__: node_class
                              for node_class in node_classes}
        return keyed_node_classes

    def generate_flow_engine(self, keyed_node_classes=None):
        flow_engine = FlowEngine()
        for node_class in keyed_node_classes.values():
            flow_engine.register_node_class(node_class=node_class)
        return flow_engine

    def generate_flow_runner(self, flow_engine=None,
                                 flow_client=None):
        flow_runner = BaseFlowRunner(flow_engine=flow_engine,
                                             flow_client=flow_client)
        return flow_runner

    def populate_flows(self):
        flow = Flow()
        root_node = self.keyed_node_classes['Node_0'](flow=flow,
                                                      status='PENDING')
        flow.add_node(node=root_node, as_root=True)
        serialization = self.flow_engine.serialize_flow(flow)
        flow_model = FlowModel.objects.create(
            serialization=json.dumps(serialization))
        self.flow_uuid = flow_model.uuid

    def test_run_flows(self):
        states = []
        states.append(self.get_state(prev_state=None))

        def tick_and_capture_state(label=None):
            print("tick_and_capture", label)
            self.flow_runner.tick()
            states.append(self.get_state(prev_state=states[-1]))

        tick_and_capture_state(label='a')
        self.assertEqual(
            states[-1]['node_classes']['Node_0']['ticks'], 1)

        tick_and_capture_state(label='b')
        self.assertEqual(
            states[-1]['node_classes']['Node_0']['ticks'], 1)
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 1)

        tick_and_capture_state(label='c')
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 2)
        self.assertFalse(
            'ticks' in states[-1]['node_classes']['Node_2'])

        self.keyed_node_classes['Node_1'].allow_completion = True
        tick_and_capture_state(label='d')
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 3)

        tick_and_capture_state(label='e')
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 3)
        self.assertEqual(
            states[-1]['node_classes']['Node_2']['ticks'], 1)

        self.keyed_node_classes['Node_2'].allow_completion = True
        tick_and_capture_state(label='f')
        self.assertEqual(
            states[-1]['node_classes']['Node_2']['ticks'], 2)
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
            'node_classes': {
                key: node_class.__dict__
                for key, node_class in self.keyed_node_classes.items()
            }
        }
        return state
