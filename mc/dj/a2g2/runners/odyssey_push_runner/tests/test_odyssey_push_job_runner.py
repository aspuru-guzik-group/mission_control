import unittest
from unittest.mock import call, Mock, patch


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
        self.runner.generate_transfer_client = Mock()
        self.runner.setup()
        self.assertEqual(self.runner.base_job_runner,
                         self.mocks['BaseJobRunner'].return_value)
        self.assertEqual(self.mocks['BaseJobRunner'].call_args,
                         call(execution_client=self.runner\
                              .generate_execution_client.return_value,
                              transfer_client=self.runner\
                              .generate_transfer_client.return_value
                             )
                        )

class GenerateExecutionClientTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['SSHClient'] = patch.object(
            odyssey_push_job_runner, 'SSHControlSocketClient')
        self.patchers['SlurmClient'] = patch.object(
            odyssey_push_job_runner, 'RemoteSlurmExecutionClient')

    def test_generates_remote_slurm_execution_client(self):
        execution_client = self.runner.generate_execution_client()
        self.assertEqual(execution_client,
                         self.mocks['SlurmClient'].return_value)
        self.assertEqual(self.mocks['SlurmClient'].call_args,
                         call(ssh_client=self.mocks['SSHClient'].return_value))
        self.assertEqual(self.mocks['SSHClient'].call_args,
                         call(user=self.runner.odyssey_user,
                              host=self.runner.odyssey_host))

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
