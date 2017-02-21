import filecmp
import json
import os
import subprocess
import tempfile
import unittest

from ..rdkit_conformer_generator import RDKitConformerGenerator
from .. import cmd

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        _dir = os.path.dirname(__file__)
        self.output_dir = tempfile.mkdtemp()
        self.expected_output_dir = os.path.join(_dir, 'expected_output')
        self.input_params = {
            'output_dir': self.output_dir,
            'forcefield_id': 'mmff',
            'smiles': ('c1ccc2c(c1)[nH]c1ccc(-c3nnc(-c4ccc5[nH]'
                       'c6ccccc6c5c4)[nH]3)cc12'),
            'output_limit': 20,
            'candidate_limit': 200,
            'prune_rms_thresh': 0.01,
            'cluster_energy_range': 5.0,
            'cluster_energy_radius': 5.0,
            'cluster_energy_auto_accept_radius': 1e-3,
            'cluster_rms_radius': 0.1,
            'xyz_filename_tpl': 'conf_{index}.xyz',
        }

    def check_output_dir(self):
        comparison = filecmp.dircmp(self.expected_output_dir, self.output_dir)
        self.assertEqual(comparison.left_list, comparison.right_list)

class GeneratorTestCase(BaseTestCase):
    def test_generates_expected_conformers(self):
        generator = RDKitConformerGenerator(**self.input_params)
        generator.generate_conformers()
        self.check_output_dir()

class CmdTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cmd_module_name = cmd.__name__

    def test_takes_params_from_cli_options(self):
        cli_options = ['--{param}="{value}"'.format(param=param, value=value)
                       for param, value in self.input_params.items()]
        cmd = "python -m %s %s" % (self.cmd_module_name, " ".join(cli_options))
        subprocess.call(cmd, shell=True)
        self.check_output_dir()

    def test_takes_params_from_file(self):
        params_file_dir = tempfile.mkdtemp(prefix='rdkit_confgen_params.')
        params_file_path = os.path.join(params_file_dir, 'params.json')
        with open(params_file_path, 'w') as f: json.dump(self.input_params, f)
        cmd = "python -m %s --params_file=%s" % (self.cmd_module_name,
                                                 params_file_path)
        subprocess.call(cmd, shell=True)
        self.check_output_dir()

if __name__ == '__main__':
    unittest.main()
