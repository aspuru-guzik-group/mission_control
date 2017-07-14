import json
import unittest

from mc.houston.houston import Houston
from mc.houston import test_utils as _houston_test_utils
from mc.utils import test_utils as _mc_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = Houston(cfg=self._generate_cfg())
        self.houston.utils.ensure_db()

    def capture(self): return _mc_test_utils.capture()
    def capture_std_streams(self, streams=None):
        return _mc_test_utils.capture_std_streams(streams=streams)

    def _generate_cfg(self):
        return _houston_test_utils.generate_test_cfg()

    def _generate_flow_spec(self, flow_spec=None):
        flow_spec = flow_spec or {
            'key': 'my_flow',
            'tasks': [
                {
                    'key': 'task_%s' % i,
                    'task_type': 'log',
                    'task_params': {'msg': 'msg_%s' % i}
                }
                for i in range(3)
            ]
        }
        return flow_spec

class TickFlowSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_spec = self._generate_flow_spec()

    def test_runs_single_tick(self):
        with self.capture_std_streams(streams=['stdout', 'stderr']) as streams:
            self.houston.call_command(
                'tick', 'flow', flow_spec=json.dumps(self.flow_spec))
        expected_flow = self.houston.utils.flow_engine.flow_spec_to_flow(
            flow_spec=self.flow_spec)
        self.houston.utils.flow_engine.tick_flow(flow=expected_flow)
        expected_flow_dict = self.houston.utils.flow_engine.flow_to_flow_dict(
            flow=expected_flow)
        output = streams['stdout'].read()
        self.assertEqual(json.loads(output), expected_flow_dict)
