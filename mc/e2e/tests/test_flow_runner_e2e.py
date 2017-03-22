import json
import unittest

from django.conf.urls import url, include
from django.test import TestCase, override_settings

from mc.mc_utils import test_utils
from missions.models import Flow as FlowModel
from mc.flow_engines.flow import Flow
from mc.flow_runners.base_flow_runner import BaseFlowRunner
from mc.flow_engines.flow_engine import FlowEngine
from mc.mc_client.mission_control_client import MissionControlClient


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('missions.urls')),
]

@unittest.skip("needs big fixup, punting right now")
@override_settings(ROOT_URLCONF=__name__)
class FlowRunnerE2ETestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client(request_client=self.client)
        self.flow_engine = self.generate_flow_engine()
        self.flow_client = MissionControlClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )
        self.flow_runner = self.generate_flow_runner(
            flow_client=self.flow_client,
            flow_engine=self.flow_engine)
        self.populate_flows()

    def generate_flow_engine(self):
        flow_engine = FlowEngine()
        return flow_engine

    def generate_flow_runner(self, flow_engine=None, flow_client=None):
        flow_runner = BaseFlowRunner(flow_engine=flow_engine,
                                     flow_client=flow_client)
        return flow_runner

    def populate_flows(self):
        flow = Flow()
        root_node = {'status': 'PENDING'}
        flow.add_node(node=root_node, as_root=True)
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
        self.assertEqual(states[-1]['node_engines']['NodeEngine_0']['ticks'], 1)

        tick_and_capture_state(label='b')
        self.assertEqual(states[-1]['node_engines']['NodeEngine_0']['ticks'], 1)
        self.assertEqual(states[-1]['node_engines']['NodeEngine_1']['ticks'], 1)

        tick_and_capture_state(label='c')
        self.assertEqual(states[-1]['node_engines']['NodeEngine_1']['ticks'], 2)
        self.assertFalse('ticks' in states[-1]['node_engines']['NodeEngine_2'])

        self.keyed_node_engines['NodeEngine_1'].allow_completion = True
        tick_and_capture_state(label='d')
        self.assertEqual(states[-1]['node_engines']['NodeEngine_1']['ticks'], 3)

        tick_and_capture_state(label='e')
        self.assertEqual(states[-1]['node_engines']['NodeEngine_1']['ticks'], 3)
        self.assertEqual(states[-1]['node_engines']['NodeEngine_2']['ticks'], 1)

        self.keyed_node_engines['NodeEngine_2'].allow_completion = True
        tick_and_capture_state(label='f')
        self.assertEqual(states[-1]['node_engines']['NodeEngine_2']['ticks'], 2)
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
            'node_engines': {
                key: node_engine.__dict__
                for key, node_engine in self.keyed_node_engines.items()
            }
        }
        return state
