import json
import os
import tempfile
import unittest
from unittest.mock import call, Mock, patch

from ..base_command import BaseCommand


class BaseCommandBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = self.generate_job()
        self.cfg = self.generate_cfg()
        self.ctx_dir = 'some_ctx_dir'
        self.output_dir = 'some_output_dir'
        self.file_args = self.generate_file_args()
        self.argv = self.generate_argv(arg_tuples=[
            *[(file_arg_name, file_arg['path']) 
              for file_arg_name, file_arg in self.file_args.items()
             ],
            ('ctx_dir', self.ctx_dir),
            ('output_dir', self.output_dir),
        ])
        self.command = self.generate_command()

    def generate_job(self): return {}

    def generate_cfg(self): return {}

    def generate_command(self):
        class BasicCommand(BaseCommand):
            def execute_job(self):
                pass
        return BasicCommand()

    def generate_file_args(self):
        file_args = {
            'job': {
                'value': self.job,
            },
            'cfg': {
                'value': self.cfg
            }
        }
        self.tmpdir = tempfile.mkdtemp()
        self.input_paths = {}
        for file_arg_name, file_arg in file_args.items():
            file_path = os.path.join(self.tmpdir, file_arg_name) +  '.json'
            with open(file_path, 'w') as f: json.dump(file_arg['value'], f)
            file_arg['path'] = file_path
        return file_args

    def generate_argv(self, arg_tuples=None):
        argv = []
        for arg_name, arg_value in arg_tuples:
            argv.append('--%s' % arg_name)
            argv.append(arg_value)
        return argv

    def _execute_command(self, command=None):
        command = command or self.command
        command.execute(argv=self.argv)

    def test_get_cfg_value_from_cfg(self):
        cfg = {'key1': 'value1'}
        cfg_value = self.command.get_cfg_value(cfg=cfg, key='key1')
        self.assertEqual(cfg_value, cfg['key1'])

    def test_get_cfg_value_from_env(self):
        cfg = {'key1': 'value1'}
        mock_environ = {'key1': 'value_from_env'}
        with patch.object(os, 'environ', new=mock_environ):
            cfg_value = self.command.get_cfg_value(cfg=cfg, key='key1')
            self.assertEqual(cfg_value, mock_environ['key1'])

    def test_calls_execute_job(self):
        self.command.execute_job = Mock()
        self._execute_command()
        self.assertEqual(self.command.execute_job.call_args,
                         self.get_expected_execute_job_call())

    def get_expected_execute_job_call(self):
        return call(job=self.job, cfg=self.cfg,
                    ctx_dir=self.ctx_dir, output_dir=self.output_dir)
