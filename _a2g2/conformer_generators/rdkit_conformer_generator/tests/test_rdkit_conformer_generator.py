import filecmp
import json
import os
import subprocess
import tarfile
import tempfile
import unittest

from ..rdkit_conformer_generator import RDKitConformerGenerator
from .. import cmd

@unittest.skipUnless(os.environ.get('MC_RUN_LONG_TESTS', False),
                     'MC_RUN_LONG_TESTS is not set')
class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
        self.output_dir = tempfile.mkdtemp()
        self.xyz_filename_tpl = 'conf_{index}.xyz'
        self.confgen_params = {
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
        }

    def check_output_dir(self):
        expected_output_dir = self.generate_expected_output_dir()
        comparison = filecmp.dircmp(self.output_dir, expected_output_dir)
        self.assertEqual(comparison.left_list, comparison.right_list)

    def generate_expected_output_dir(self):
        output_dir = tempfile.mkdtemp()
        xyzs_tgz_path = os.path.join(self.fixtures_dir, 'xyzs.tar.gz')
        conformers_outdir = os.path.join(
            output_dir, RDKitConformerGenerator.CONFORMERS_OUTDIR_NAME)
        tarfile.open(xyzs_tgz_path).extractall(conformers_outdir)
        params_outfile_path = os.path.join(
            output_dir, RDKitConformerGenerator.PARAMS_OUTFILE_NAME)
        json.dump(self.confgen_params, open(params_outfile_path, 'w'))
        return output_dir

class GeneratorTestCase(BaseTestCase):
    def test_generates_expected_conformers(self):
        generator = RDKitConformerGenerator(
            output_dir=self.output_dir,
            xyz_filename_tpl=self.xyz_filename_tpl,
            **self.confgen_params)
        generator.generate_conformers()
        self.check_output_dir()

class CmdTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cmd_module_name = cmd.__name__
        self.cli_params = {
            'output_dir': self.output_dir,
            'xyz_filename_tpl': self.xyz_filename_tpl,
            **self.confgen_params,
        }

    def test_takes_params_from_cli_options(self):
        cli_options = [
            '--{param}="{value}"'.format(param=param, value=value)
            for param, value in self.cli_params.items()
        ]
        cmd = "python -m %s %s" % (self.cmd_module_name, " ".join(cli_options))
        subprocess.call(cmd, shell=True)
        self.check_output_dir()

    def test_takes_params_from_file(self):
        params_file_dir = tempfile.mkdtemp(prefix='rdkit_confgen_params.')
        params_file_path = os.path.join(params_file_dir, 'params.json')
        with open(params_file_path, 'w') as f: json.dump(self.cli_params, f)
        cmd = "python -m %s --params_file=%s" % (self.cmd_module_name,
                                                 params_file_path)
        subprocess.call(cmd, shell=True)
        self.check_output_dir()

if __name__ == '__main__':
    unittest.main()
