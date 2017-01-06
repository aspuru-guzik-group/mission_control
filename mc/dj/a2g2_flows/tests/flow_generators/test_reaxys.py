import collections
import unittest
from unittest.mock import MagicMock

from flow_engines.flow_engine import FlowEngine

from ...flow_generators.reaxys import ReaxysFlowGenerator
from ...flow_generators.reaxys import RunAndLoadFlowGenerator


class ReaxysFlowTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx = collections.defaultdict(MagicMock)
        self.engine = FlowEngine(ctx=self.ctx)
        self.flow_spec = {
            'type': 'Reaxys',
            'smiles': 'some smiles',
            'confgen_params': 'some confgen params',
        }
        self.flow = ReaxysFlowGenerator.generate_flow(flow_spec=self.flow_spec)

    def test_run_flow(self):
        # Helpers
        self.tick_counter = 0
        def _tick():
            self.tick_counter += 1
            self.engine.tick_flow(flow=self.flow, ctx=self.ctx)
        tasks = self.flow.tasks

        # t = 0

        # t = 1
        _tick()
        self.assertEqual(tasks['confgen'].status, 'RUNNING')
        self.assertEqual(tasks['confgen'].input, {
            'flow_spec': {
                'flow_type': RunAndLoadFlowGenerator.flow_type,
                'run_spec': {
                    'job_type': 'confgen',
                    'confgen': {
                        'smiles': self.flow_spec['smiles'],
                        'params': self.flow_spec['confgen_params'],
                    }
                },
                'load_spec': {
                    'job_type': 'confgen:load',
                }
            }
        })

        confgen_flow = {
            'status': 'COMPLETED',
            'output': {'dir': 'output dir'}
        }
        self.ctx['get_flow'].return_value = confgen_flow

        # t = 2
        _tick()
        self.assertEqual(tasks['confgen'].status, 'COMPLETED')
        self.assertEqual(self.flow.output, confgen_flow['output'])

    def get_state(self, prev_state=None):
        state = {}
        return state

if __name__ == '__main__':
    unittest.main()
