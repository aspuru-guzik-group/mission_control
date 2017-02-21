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
                'type': 'job_type'
            }
        }
        self.cfg = defaultdict(MagicMock)

class ExecuteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a2g2_job_engine.get_job_module = MagicMock()

    def test_dispatches_to_get_job_module_result(self):
        result = self.a2g2_job_engine.execute_job(job=self.job, cfg=self.cfg)
        self.assertEqual(self.a2g2_job_engine.get_job_module.call_args,
                         call(job=self.job, cfg=self.cfg))
        expected_module = self.a2g2_job_engine.get_job_module.return_value
        self.assertEqual(expected_module.execute_job.call_args,
                         call(job=self.job, cfg=self.cfg))
        self.assertEqual(result, expected_module.execute_job.return_value)

class GetJobModuleCase(BaseTestCase):
    @patch.object(a2g2_job_engine, 'importlib')
    def test_imports_expected_module(self, patched_importlib):
        result = self.a2g2_job_engine.get_job_module(job=self.job, cfg=self.cfg)
        expected_module_name = '{job_type}_job_module'.format(
            job_type=self.job['job_spec']['type'])
        self.assertEqual(
            patched_importlib.import_module.call_args,
            call(expected_module_name, package=a2g2_job_engine.__name__))
        self.assertEqual(result, patched_importlib.import_module.return_value)

class ExecuteJobCommandTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cmd = a2g2_job_engine.ExecuteJobCommand()

    def generate_argv(self, arg_files=None):
        argv = []
        for arg_name, arg_file in arg_files.items():
            argv.append('--%s' % arg_name)
            argv.append(arg_file)
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

    def test_calls_execute_job_w_parsed_args(self):
        mock_generate_a2g2_job_engine = MagicMock()
        self.cmd.generate_a2g2_job_engine = mock_generate_a2g2_job_engine
        arg_files = self.generate_arg_files()
        self.cmd.execute(argv=self.generate_argv(arg_files=arg_files))
        self.assertEqual(
            mock_generate_a2g2_job_engine.return_value.execute_job.call_args,
            call(**{arg_name: json.load(open(arg_file_path))
                    for arg_name, arg_file_path in arg_files.items()})
        )

