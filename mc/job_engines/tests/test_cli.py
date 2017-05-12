import json
import os
import shutil
import tempfile
import unittest
from unittest.mock import call, MagicMock

from .. import cli


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.cmd = cli.JobEngineCommand()
        self.job_engine = MagicMock()
        self.cmd.get_job_engine = MagicMock(return_value=self.job_engine)

    def generate_argv(self, job_command=None, arg_tuples=None):
        argv = [job_command]
        for arg_name, arg_value in arg_tuples:
            argv.append('--%s' % arg_name)
            argv.append(arg_value)
        return argv

class BuildSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.tmp_dir = tempfile.mkdtemp()
        self.job = {'some': 'job'}
        self.cfg = {'some': 'cfg'}
        self.output_dir = os.path.join(self.tmp_dir, 'output')
        self.arg_files = self.generate_arg_files(arg_names=['job', 'cfg'])
        self.arg_tuples = [*self.arg_files.items(),
                           ('output_dir', self.output_dir)]
        self.job_command = 'build_submission'
        self.argv = self.generate_argv(job_command=self.job_command,
                                       arg_tuples=self.arg_tuples)
    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def generate_arg_files(self, arg_names=None):
        arg_files = {}
        for arg_name in arg_names:
            arg_file_path = os.path.join(self.tmp_dir, arg_name + '.json')
            with open(arg_file_path, 'w') as f:
                json.dump(getattr(self, arg_name), f)
            arg_files[arg_name] = arg_file_path
        return arg_files

    def test_calls_execute_command_w_parsed_args(self):
        self.cmd.execute(argv=self.argv)
        expected_command_kwargs = { 
            'command': self.job_command,
            **{arg_name: json.load(open(arg_file_path))
               for arg_name, arg_file_path in self.arg_files.items()},
            'output_dir': self.output_dir,
        }
        self.assertEqual(self.job_engine.execute_command.call_args,
                         call(**expected_command_kwargs))

class RunSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_dir = 'some_submission_dir'
        self.arg_tuples = [('submission_dir', self.submission_dir)]
        self.job_command = 'run_submission'
        self.argv = self.generate_argv(job_command=self.job_command,
                                       arg_tuples=self.arg_tuples)

    def test_calls_execute_command_w_parsed_args(self):
        self.cmd.execute(argv=self.argv)
        expected_command_kwargs = {
            'command': self.job_command,
            'submission_dir': self.submission_dir
        }
        self.assertEqual(self.job_engine.execute_command.call_args,
                         call(**expected_command_kwargs))
