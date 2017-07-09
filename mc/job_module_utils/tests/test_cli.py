import unittest
from unittest.mock import call, MagicMock

from .. import cli


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.cmd = cli.JobModuleCliCommand()
        self.cmd._get_command_dispatcher = MagicMock()

    def generate_argv(self, job_command=None, arg_tuples=None):
        argv = ['cmd_exe', job_command]
        for arg_name, arg_value in arg_tuples:
            argv.append('--%s' % arg_name)
            argv.append(arg_value)
        return argv

class RunJobdirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_command = 'run_jobdir'
        self.jobdir = 'some_jobdir'
        self.arg_tuples = [('jobdir', self.jobdir)]
        self.argv = self.generate_argv(job_command=self.job_command,
                                       arg_tuples=self.arg_tuples)

    def test_calls_command_w_parsed_args(self):
        self.cmd.run_from_argv(self.argv)
        expected_command_kwargs = {'jobdir': self.jobdir}
        self.assertEqual(
            self.cmd._get_command_dispatcher.return_value.run_jobdir.call_args,
            call(**expected_command_kwargs)
        )
