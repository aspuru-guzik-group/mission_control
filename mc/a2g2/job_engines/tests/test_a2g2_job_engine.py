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

class ExecuteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a2g2_job_engine.get_job_engine = MagicMock()

    def test_dispatches_to_get_job_engine_result(self):
        result = self.a2g2_job_engine.execute_job(
            job=self.job,
            cfg=self.cfg,
            output_dir=self.output_dir,
            ctx_dir=self.ctx_dir
        )
        self.assertEqual(self.a2g2_job_engine.get_job_engine.call_args,
                         call(job=self.job, cfg=self.cfg))
        expected_engine = self.a2g2_job_engine.get_job_engine.return_value
        self.assertEqual(
            expected_engine.execute_job.call_args,
            call(job=self.job, cfg=self.cfg, output_dir=self.output_dir,
                 ctx_dir=self.ctx_dir)
        )
        self.assertEqual(result, expected_engine.execute_job.return_value)

class GetJobModuleCase(BaseTestCase):
    @patch.object(a2g2_job_engine, 'importlib')
    def test_imports_expected_module(self, patched_importlib):
        result = self.a2g2_job_engine.get_job_module(job=self.job, cfg=self.cfg)
        # @TODO: brittle test, because of packages, make more robust later.
        expected_module_name = '{package}.{module}'.format(
            package=a2g2_job_engine.DEFAULT_JOB_MODULE_PKGS[0],
            module=self.job['job_spec']['module']
        )
        self.assertEqual(patched_importlib.import_module.call_args,
                         call(expected_module_name))
        self.assertEqual(result, patched_importlib.import_module.return_value)

class ExecuteJobCommandTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cmd = a2g2_job_engine.ExecuteJobCommand()
        self.arg_files = self.generate_arg_files()
        self.arg_tuples = [
            *self.arg_files.items(),
            ('output_dir', self.output_dir),
            ('ctx_dir', self.ctx_dir)
        ]
        self.argv = self.generate_argv(arg_tuples=self.arg_tuples)

    def generate_argv(self, arg_tuples=None):
        argv = []
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

    def test_calls_execute_job_w_parsed_args(self):
        mock_generate_job_engine = MagicMock()
        self.cmd.generate_job_engine = mock_generate_job_engine
        self.cmd.execute(argv=self.argv)
        self.assertEqual(
            mock_generate_job_engine.return_value.execute_job.call_args,
            call(
                **{
                    **{
                        arg_name: json.load(open(arg_file_path))
                        for arg_name, arg_file_path in self.arg_files.items()
                    },
                    'output_dir': self.output_dir,
                    'ctx_dir': self.ctx_dir
                }
            )
        )
