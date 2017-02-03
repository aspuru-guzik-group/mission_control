import os
import tempfile
import unittest
from unittest.mock import call, Mock

from rdkit import AllChem, Chem

from ..confgen_load_job_engine import ConfgenLoadJobEngine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = ConfgenLoadJobEngine()
        self.conformers = self.generate_conformers()

    def generate_conformers(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('Cc1ccccc1'))
        conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=10)
        conformers = mol.GetConformers()

class RunJobDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_dir = self.generate_job_dir()
        self.job_spec = self.generate_job_spec()

    def generate_job_dir(self):
        job_dir = tempfile.mkdtemp()
        raw_data_dir = self.generate_raw_data_dir(parent_dir=job_dir)
        for i, conformer in enumerate(self.conformers):
            fmt = 'pdb'
            conformer_file_name = 'conformer_%s.%s' % (i, fmt)
            conformer_file_path = os.path.join(raw_data_dir,
                                               conformer_file_name)
            self.output_conformer(conformer=conformer,
                                  format='pdb',
                                  output_file_path=conformer_file_path)
        return job_dir

    def generate_raw_data_dir(self, parent_dir=None):
        raw_data_dir = os.path.join(parent_dir, 'raw_data')
        os.makedirs(raw_data_dir)

    def test_uploads_parse_results(self):
        self.engine.parse_inputs = Mock()
        self.engine.upload_objects = Mock()
        self.engine.execute_job(job_dir=self.job_dir)
        self.assertEqual(
            self.engine.parse_inputs.call_args,
            call(inputs_dir=os.path.join(self.job_dir, 'raw_data'))
        )
        self.assertEqual(self.engine.upload_objects.call_args,
                         call(objects=self.engine.parse_inputs.return_value))

class ParseInputsTestCase(BaseTestCase):
    def test_generates_objects_from_inputs(self):
        self.fail()

class UploadTestCase(BaseTestCase):
    def test_makes_expected_uploads(self):
        self.fail()

    def test_outputs_upload_results(self):
        self.fail()
