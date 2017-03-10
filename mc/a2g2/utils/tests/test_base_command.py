import json
import os
import tempfile
import unittest
from unittest.mock import patch

from ..base_command import BaseCommand


class BaseCommandBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.cfg = self.generate_cfg()
        self.tmpdir = tempfile.mkdtemp()
        self.file_args = self.generate_file_args(target_dir=self.tmpdir)
        self.argv = self.generate_argv(arg_tuples=[
            *[(file_arg_name, file_arg['path']) 
              for file_arg_name, file_arg in self.file_args.items()
             ],
        ])
        self.command = self.generate_command()

    def generate_cfg(self): return {}

    def generate_command(self):
        class BasicCommand(BaseCommand): pass
        return BasicCommand()

    def generate_file_args(self, target_dir=None):
        file_args = {
            'cfg': {
                'value': self.cfg
            }
        }
        self.input_paths = {}
        for file_arg_name, file_arg in file_args.items():
            file_arg['path'] = self.generate_json_file_for_arg(
                target_dir=target_dir,
                arg_name=file_arg_name,
                arg_value=file_arg['value']
            )
        return file_args

    def generate_json_file_for_arg(self, target_dir=None, arg_name=None,
                                   arg_value=None):
        file_path = os.path.join(target_dir, arg_name) +  '.json'
        with open(file_path, 'w') as f: json.dump(arg_value, f)
        return file_path

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
