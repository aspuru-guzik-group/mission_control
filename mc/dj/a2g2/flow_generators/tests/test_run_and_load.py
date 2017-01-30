from collections import defaultdict
import logging
import unittest
from unittest.mock import call, MagicMock

from flow_engines.flow_engine import FlowEngine

from ..run_and_load import RunAndLoadFlowGenerator


class RunAndLoadFlowTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx = self.generate_ctx()
        self.engine = FlowEngine(ctx=self.ctx)
        self.flow_spec = {
            'run_spec': defaultdict(MagicMock),
            'load_spec': defaultdict(MagicMock),
        }
        self.flow = RunAndLoadFlowGenerator.generate_flow(
            flow_spec=self.flow_spec)

    def generate_ctx(self):
        return {
            'create_job': MagicMock(),
            'get_job': MagicMock(),
        }

    def test_flow(self):
        # Helpers
        self.tick_counter = 0
        def _tick():
            logging.info("tick_counter: %s" % self.tick_counter)
            self.tick_counter += 1
            self.engine.tick_flow(flow=self.flow, ctx=self.ctx)
        tasks = self.flow.tasks

        # t = 0
        self.assertEqual(tasks['run'].status, 'PENDING')

        # t = 1
        _tick()
        self.assertEqual(tasks['run'].input['job_spec'],
                         self.flow_spec['run_spec'])
        self.assertEqual(tasks['run'].status, 'RUNNING')
        self.assertEqual(
            self.ctx['create_job'].call_args,
            call(job_kwargs={'spec': tasks['run'].input['job_spec']}))

        mock_run_job = {'status': 'COMPLETED',
                        'data': {'output': {'raw_dir': 'some dir'}}}
        self.ctx['get_job'].return_value = mock_run_job

        # t = 2
        _tick()
        self.assertEqual(tasks['run'].status, 'COMPLETED')
        self.assertEqual(tasks['run'].output, mock_run_job['data']['output'])
        self.assertEqual(tasks['load_prep'].status, 'PENDING')

        # t = 3
        _tick()
        self.assertEqual(tasks['load_prep'].status, 'COMPLETED')
        self.assertEqual(tasks['load'].status, 'PENDING')

        # t = 4
        _tick()
        self.assertEqual(tasks['load'].status, 'RUNNING')
        self.assertEqual(
            self.ctx['create_job'].call_args,
            call(job_kwargs={ 
                'spec': {
                    **self.flow_spec['load_spec'],
                    'raw_dir': mock_run_job['data']['output']['raw_dir']
                }
            })
        )

        mock_load_job = {'status': 'COMPLETED',
                         'data': {'output': 'some output'}}
        self.ctx['get_job'].return_value = mock_load_job

        # t = 5
        _tick()
        self.assertEqual(self.flow.status, 'COMPLETED')
        self.assertEqual(self.flow.output, mock_load_job['data']['output'])

if __name__ == '__main__':
    unittest.main()
