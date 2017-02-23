from collections import defaultdict
import os
import unittest
from unittest.mock import call, MagicMock, Mock, patch

from .. import odyssey_push_job_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = self.start_patchers()
        self.runner = self.generate_runner()

    def decorate_patchers(self): pass

    def start_patchers(self):
        mocks = {key: patcher.start() for key, patcher in self.patchers.items()}
        return mocks

    def tearDown(self):
        self.stop_patchers()

    def stop_patchers(self):
        for patcher in self.patchers.values(): patcher.stop()

    def generate_runner(self, *args, **kwargs):
        default_kwargs = {'setup': False, 'odyssey_user': Mock(),
                          'odyssey_host': Mock()}
        return odyssey_push_job_runner.OdysseyPushJobRunner(
            *args, **{**default_kwargs, **kwargs})

class SetupTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['BaseJobRunner'] = patch.object(odyssey_push_job_runner,
                                                      'BaseJobRunner')

    def test_generates_base_job_runner(self):
        self.runner.generate_execution_client = Mock()
        self.runner.setup()
        self.assertEqual(self.runner.base_job_runner,
                         self.mocks['BaseJobRunner'].return_value)
        self.assertEqual(
            self.mocks['BaseJobRunner'].call_args,
            call(execution_client=\
                 self.runner.generate_execution_client.return_value,
                 get_job_execution_result=\
                 self.runner.get_job_execution_result)
        )

class GenerateExecutionClientTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['ExecutionClient'] = patch.object(
            odyssey_push_job_runner, 'ExecutionClient')

    def test_generates_remote_slurm_execution_client(self):
        ssh_client = Mock()
        execution_client = self.runner.generate_execution_client(
            ssh_client=ssh_client)
        self.assertEqual(execution_client,
                         self.mocks['ExecutionClient'].return_value)
        self.assertEqual(self.mocks['ExecutionClient'].call_args,
                         call(ssh_client=ssh_client))

class GetJobExecutionResultTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_state = defaultdict(MagicMock)

    @patch('os.path.exists')
    def test_returns_success_if_completed_checkpoint_file_exists(
        self, mock_exists):
        mock_exists.return_value = True
        result = self.runner.get_job_execution_result(
            job_state=self.job_state)
        expected_checkpoint_path = os.path.join(
            self.job_state['completed_dir'],
            self.job_state['submission']['checkpoint_files']['completed']
        )
        self.assertEqual(mock_exists.call_args,
                         call(expected_checkpoint_path))
        self.assertEqual(result, {'result': 'COMPLETED'})

    @patch.object(odyssey_push_job_runner, 'open')
    @patch('os.path.exists')
    def test_returns_error_if_job_failed(self, mock_exists, mock_open):
        mock_exists.return_value = False
        result = self.runner.get_job_execution_result(job_state=self.job_state)
        expected_checkpoint_path = os.path.join(
            self.job_state['completed_dir'],
            self.job_state['submission']['checkpoint_files']['failed']
        )
        self.assertEqual(mock_open.call_args, call(expected_checkpoint_path))
        self.assertEqual(
            result,
            {
                'result': 'FAILED',
                'error': mock_open.return_value.read.return_value
            }
        )

class RunTestCase(BaseTestCase):
    def test_wraps_base_job_runner(self):
        self.runner.base_job_runner = Mock()
        self.runner.run()
        self.assertEqual(self.runner.base_job_runner.run.call_count, 1)

class TickTestCase(BaseTestCase):
    def test_wraps_base_job_runner(self):
        self.runner.base_job_runner = Mock()
        self.runner.tick()
        self.assertEqual(self.runner.base_job_runner.tick.call_count, 1)

if __name__ == '__main__':
    unittest.main()
