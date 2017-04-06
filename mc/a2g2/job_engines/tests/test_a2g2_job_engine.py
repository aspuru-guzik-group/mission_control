from collections import defaultdict
import json
import os
import tempfile
import unittest
from unittest.mock import call, patch, MagicMock

from .. import a2g2_job_engine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.a2g2_job_engine = a2g2_job_engine.A2G2JobEngine()
        self.job = {
            'job_spec': {
                'module': 'job_module'
            }
        }
        self.cfg = defaultdict(MagicMock)
        self.output_dir = 'some_output_dir'
        self.ctx_dir = 'some_ctx_dir'

class BuildJobSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.output_dir = MagicMock()
        self.a2g2_job_engine.get_job_module = MagicMock()
        self.expected_module = self.a2g2_job_engine.get_job_module.return_value

    def test_dispatches_to_job_module(self):
        self.a2g2_job_engine.build_job_submission(
            job=self.job, cfg=self.cfg, output_dir=self.output_dir)
        self.assertEqual(self.a2g2_job_engine.get_job_module.call_args,
                         call(job=self.job, cfg=self.cfg))
        self.assertEqual(
            self.expected_module.build_job_submission.call_args,
            call(job=self.job, cfg=self.cfg, job_dir=self.output_dir))

class ExecuteCommandTestCase(BaseTestCase):
    def test_build_job_submission_command(self):
        self.assert_handles_command(command='build_job_submission')

    def assert_handles_command(self, command=None):
        mock_handler = MagicMock()
        setattr(self.a2g2_job_engine, command, mock_handler)
        kwargs = {'job': self.job, 'cfg': self.cfg, 'extra': 'extra'}
        result = self.a2g2_job_engine.execute_command(command=command, **kwargs)
        self.assertEqual(mock_handler.call_args, call(**kwargs))
        self.assertEqual(result, mock_handler.return_value)

    def test_run_job_submission_command(self):
        self.assert_handles_command(command='run_job_submission')

class GetJobModuleCase(BaseTestCase):
    @patch.object(a2g2_job_engine, 'importlib')
    def test_imports_expected_module(self, patched_importlib):
        self.a2g2_job_engine.get_job_module_name = MagicMock(
            return_value='some_module')
        result = self.a2g2_job_engine.get_job_module(job=self.job, cfg=self.cfg)
        # @TODO: brittle test, because of packages, make more robust later.
        expected_module_name = '{package}.{module}'.format(
            package=a2g2_job_engine.DEFAULT_JOB_MODULE_PKGS[0],
            module=self.a2g2_job_engine.get_job_module_name.return_value,
        )
        self.assertEqual(patched_importlib.import_module.call_args,
                         call(expected_module_name))
        self.assertEqual(result, patched_importlib.import_module.return_value)

class JobEngineCommandTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cmd = a2g2_job_engine.JobEngineCommand()
        self.arg_files = self.generate_arg_files()
        self.arg_tuples = [
            *self.arg_files.items(),
            ('output_dir', self.output_dir),
        ]
        self.job_command = 'build_job_submission'
        self.argv = self.generate_argv(
            job_command=self.job_command,
            arg_tuples=self.arg_tuples
        )

    def generate_argv(self, job_command=None, arg_tuples=None):
        argv = [job_command]
        for arg_name, arg_value in arg_tuples:
            argv.append('--%s' % arg_name)
            argv.append(arg_value)
        return argv

    def generate_arg_files(self):
        arg_files = {}
        tmp_dir = tempfile.mkdtemp()
        for arg_name in ['job', 'cfg']:
            arg_file_path = os.path.join(tmp_dir, arg_name + '.json')
            with open(arg_file_path, 'w') as f:
                json.dump(getattr(self, arg_name), f)
            arg_files[arg_name] = arg_file_path
        return arg_files

    def test_calls_execute_command_w_parsed_args(self):
        mock_generate_job_engine = MagicMock()
        self.cmd.generate_job_engine = mock_generate_job_engine
        self.cmd.execute(argv=self.argv)
        expected_command_kwargs = { 
            'command': self.job_command,
            **{arg_name: json.load(open(arg_file_path))
               for arg_name, arg_file_path in self.arg_files.items()},
            'output_dir': self.output_dir,
        }
        self.assertEqual(
            mock_generate_job_engine.return_value.execute_command.call_args,
            call(**expected_command_kwargs)
        )
