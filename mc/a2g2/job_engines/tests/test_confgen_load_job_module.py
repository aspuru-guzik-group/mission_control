import json
import os
import tempfile
import unittest
from unittest.mock import call, patch, Mock

import pybel
from rdkit import Chem
from rdkit.Chem import AllChem

from .. import confgen_load_job_module


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
        self.mol = self.generate_mol_w_conformers()
        self.job_dir = self.generate_job_dir()

    def generate_mol_w_conformers(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('Cc1ccccc1'))
        AllChem.EmbedMultipleConfs(mol, numConfs=3)
        return mol

    def generate_job_dir(self):
        job_dir = tempfile.mkdtemp()
        self.generate_raw_data_dir(parent_dir=job_dir)
        return job_dir

    def generate_raw_data_dir(self, parent_dir=None):
        raw_data_dir = os.path.join(parent_dir, 'raw_data')
        os.makedirs(raw_data_dir)
        conformer_file_path = os.path.join(raw_data_dir, 'conformers.sdf')
        with open(conformer_file_path, 'w') as conformer_file:
            writer = Chem.SDWriter(conformer_file)
            for conformer in self.mol.GetConformers():
                writer.write(self.mol, confId=conformer.GetId())
            writer.close()

    def test_generates_chemthings_from_inputs(self):
        chemthings = self.engine.parse_job_dir(self.job_dir,
                                               includes=['raw_data/*'])
        expected_chemthings = self.generate_chemthings_for_mol(mol=self.mol)
        self.assertEqual(self.get_sorted_chemthings(chemthings),
                         self.get_sorted_chemthings(expected_chemthings))

    def generate_chemthings_for_mol(self, mol=None):
        chemthings = [
            self.generate_chemthing_for_conformer(conformer,
                                                  props={'mol_idx': i})
            for i, conformer in enumerate(mol.GetConformers())
        ]
        return chemthings

    def generate_chemthing_for_conformer(self, conformer=None, props=None):
        chemthing = {
            'cml': self.generate_cml_for_conformer(conformer),
            'props': props,
        }
        return chemthing

    def generate_cml_for_conformer(self, conformer=None):
        serialized_conformer = Chem.MolToMolBlock(conformer.GetOwningMol(),
                                                  confId=conformer.GetId())
        pb_mol = pybel.readstring('sdf', serialized_conformer)
        return pb_mol.write('cml')

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

class CmdTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.file_args = self.generate_file_args()
        self.argv = self.generate_argv(arg_tuples=[
            *[(file_arg_name, file_arg['path']) 
              for file_arg_name, file_arg in self.file_args.items()
             ],
            ('ctx_dir', self.ctx_dir)
        ])
        self.command = confgen_load_job_module.Command()

    def generate_file_args(self):
        file_args = {
            'job': {
                'value': {}
            },
            'cfg': {
                'value': {
                    'A2G2_CLIENT_CFG_JSON': json.dumps({'base_url': 'mock_url'})
                }
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

    def _execute_command(self):
        self.command.execute(argv=self.argv)

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

    def test_generates_a2g2_client(self):
        with patch.object(confgen_load_job_module, 'A2G2_Client') as MockClient:
            with patch.object(self.command, 'execute_job'):
                self._execute_command()
                self.assertEqual(
                    MockClient.call_args,
                    call(**json.loads(
                        self.file_args['cfg']['value']['A2G2_CLIENT_CFG_JSON']))
                )

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
                    call(job=self.file_args['job']['value'],
                         ctx_dir=self.ctx_dir)
                )

