import json
import os
import tempfile
import unittest
from unittest.mock import call, patch, Mock

from .. import confgen_load_job_module
from .. import chemthing_utils
from .test_base_command import BaseCommandBaseTestCase


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx_dir = 'some_ctx_dir'
        self.engine = confgen_load_job_module.ConfgenLoadJobEngine()

class ExecuteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.engine.parse_job_dir = Mock()
        self.engine.upload_chemthings = Mock()
        self.job = {
            'data': {
                'input': {
                    'dir_to_parse': 'dir_to_parse'
                }
            }
        }
        self.engine.execute_job(job=self.job, ctx_dir=self.ctx_dir)


    def test_calls_parse_on_dir_to_parse(self):
        expected_dir_to_parse = os.path.join(
            self.ctx_dir, self.job['data']['input']['dir_to_parse'])
        self.assertEqual(
            self.engine.parse_job_dir.call_args,
            call(dir_to_parse=expected_dir_to_parse)
        )

    def test_uploads_parse_results(self):
        self.assertEqual(self.engine.upload_chemthings.call_args,
                         call(self.engine.parse_job_dir.return_value))

    def test_sets_job_output_to_upload_results(self):
        self.assertEqual(self.job['data']['output'],
                         self.engine.upload_chemthings.return_value)

class ParseJobDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mol = chemthing_utils.generate_mols_w_conformers(
            smiles_list=['Cc1ccccc1'])[0]
        self.job_dir = self.generate_job_dir()

    def generate_job_dir(self):
        job_dir = tempfile.mkdtemp()
        self.generate_raw_data_dir(parent_dir=job_dir)
        return job_dir

    def generate_raw_data_dir(self, parent_dir=None):
        raw_data_dir = os.path.join(parent_dir, 'raw_data')
        os.makedirs(raw_data_dir)
        chemthing_utils.write_conformers_to_sd_file(
            conformers=self.mol.GetConformers(),
            sd_file_path=os.path.join(raw_data_dir, 'conformers.sdf')
        )

    def test_generates_chemthings_from_inputs(self):
        chemthings = self.engine.parse_job_dir(self.job_dir,
                                               includes=['raw_data/*'])
        expected_chemthings = self.generate_chemthings_for_mol(mol=self.mol)
        self.assertEqual(self.get_sorted_chemthings(chemthings),
                         self.get_sorted_chemthings(expected_chemthings))

    def generate_chemthings_for_mol(self, mol=None):
        chemthings = [
            chemthing_utils.generate_chemthing_for_conformer(
                conformer=conformer, props={'mol_idx': i})
            for i, conformer in enumerate(mol.GetConformers())
        ]
        return chemthings

    def get_sorted_chemthings(self, chemthings):
        return sorted(chemthings, key=lambda cthing: cthing['props']['mol_idx'])

class UploadChemthingsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.chemthings = [Mock() for i in range(3)]
        self.engine.a2g2_client = Mock()
        self.result = self.engine.upload_chemthings(chemthings=self.chemthings)

    def test_creates_chemthings_via_client(self):
        expected_call_args_list = [call(chemthing=chemthing)
                                   for chemthing in self.chemthings]
        self.assertEqual(
            self.engine.a2g2_client.create_chemthing.call_args_list,
            expected_call_args_list
        )

    def test_outputs_upload_results(self):
        expected_result = [self.engine.a2g2_client.create_chemthing.return_value
                           for chemthing in self.chemthings]
        self.assertEqual(self.result, expected_result)


class CommandTestCase(BaseCommandBaseTestCase):
    def generate_command(self): return confgen_load_job_module.Command()

    def generate_cfg(self):
        return {'A2G2_CLIENT_CFG_JSON': '{}'}

    def test_calls_execute_job(self):
        with patch.object(self.command, 'generate_a2g2_client') \
                as mock_gen_client:
            with patch.object(confgen_load_job_module, 'ConfgenLoadJobEngine') \
                    as MockEngine:
                self._execute_command()
                self.assertEqual(
                    MockEngine.call_args,
                    call(a2g2_client=mock_gen_client.return_value)
                )
                self.assertEqual(
                    MockEngine.return_value.execute_job.call_args,
                    call(job=self.job, ctx_dir=self.ctx_dir)
                )

    def test_generates_a2g2_client(self):
        with patch.object(confgen_load_job_module, 'A2G2_Client') as MockClient:
            self.command.generate_a2g2_client(cfg=self.cfg)
            self.assertEqual(
                MockClient.call_args,
                call(**json.loads(self.cfg['A2G2_CLIENT_CFG_JSON']))
            )
