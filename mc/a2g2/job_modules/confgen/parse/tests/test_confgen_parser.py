import glob
import json
import os
import tempfile
import unittest

from mc.mc_utils import test_utils as mc_test_utils
from mc.a2g2.utils import rdkit_utils

from ..confgen_parser import ConfgenParser


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.parser = ConfgenParser()

class ParseConfgenDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.completed_confgen_dir  = self.generate_completed_confgen_dir()
        self.output_dir = tempfile.mkdtemp()
        self.expected_output_dir = self.generate_expected_output_dir(
            completed_confgen_dir=self.completed_confgen_dir)

    def generate_completed_confgen_dir(self):
        output_dir = tempfile.mkdtemp()
        self.generate_job_file(parent_dir=output_dir)
        self.generate_conformer_dir(parent_dir=output_dir)
        return output_dir

    def generate_job_file(self, parent_dir=None):
        job = {
            'job_spec': {
                'command': 'generate_conformers',
                'parameters': {
                    'kwarg1': 'value1'
                },
                'precursors': ['some_precursor_uri'],
            }
        }
        with open(os.path.join(parent_dir, 'job.json'), 'w') as job_file:
            json.dump(job, job_file)
        return job_file

    def generate_conformer_dir(self, parent_dir=None):
        conformer_dir = os.path.join(parent_dir, 'conformers')
        os.mkdir(conformer_dir)
        mol = rdkit_utils.generate_mols_w_conformers(
            smiles_list=['Cc1ccccc1'])[0]
        rdkit_utils.write_conformers_to_xyz_files(
            conformers=mol.GetConformers(),
            parent_dir=conformer_dir
        )
        return conformer_dir

    def generate_expected_output_dir(self, completed_confgen_dir=None):
        output_dir = tempfile.mkdtemp()
        calc_chemthing = self.get_calc_chemthing(
            completed_confgen_dir=completed_confgen_dir)
        conformer_chemthings = self.get_conformer_chemthings(
            completed_confgen_dir=completed_confgen_dir,
            calc_chemthing=calc_chemthing)
        self.write_chemthings_bulk_file(
            chemthings=[calc_chemthing, *conformer_chemthings],
            target_path=os.path.join(output_dir, 'chemthings.bulk')
        )
        return output_dir

    def get_calc_chemthing(self, completed_confgen_dir=None):
        job = self.parse_job_file(completed_confgen_dir=completed_confgen_dir)
        calc_chemthing = {
            'uri': 'a2g2:calculation:some_uri',
            'types': ['a2g2:calculation', 'a2g2:calculation:confgen'],
            'props': {
                'a2g2:confgen:parameters': job['job_spec']['parameters']
            },
            'precursors': [job['job_spec']['precursors']]
        }
        return calc_chemthing

    def parse_job_file(self, completed_confgen_dir=None):
        job_file_path = os.path.join(completed_confgen_dir, 'job.json')
        with open(job_file_path, 'r') as job_file: job = json.load(job_file)
        return job

    def get_conformer_chemthings(self, completed_confgen_dir=None,
                                 calc_chemthing=None):
        conformer_chemthings = []
        xyz_paths = glob.glob(completed_confgen_dir + '/conformers/*')
        for i, xyz_path in enumerate(xyz_paths):
            with open(xyz_path, 'r') as xyz_file: xyz = xyz_file.read()
            conformer_chemthing = {
                'uri': 'a2g2:molecule3d:uri_%s' % i,
                'types': ['a2g2:molecule', 'a2g2:molecule3d'],
                'props': {
                    'a2g2:xyz': xyz
                },
                'precursors': calc_chemthing['uri'],
            }
            conformer_chemthings.append(conformer_chemthing)
        return conformer_chemthings

    def write_chemthings_bulk_file(self, chemthings=None, target_path=None):
        with open(target_path, 'w') as bulk_file:
            for chemthing in chemthings:
                bulk_file.write(json.dumps(chemthing))

    def test_generates_expected_output_dir(self):
        self.parser.parse_completed_confgen_dir(
            completed_confgen_dir=self.completed_confgen_dir,
            output_dir=self.output_dir
        )
        mc_test_utils.assert_dirs_equal(
            test=self,
            left=self.output_dir,
            right=self.expected_output_dir
        )
