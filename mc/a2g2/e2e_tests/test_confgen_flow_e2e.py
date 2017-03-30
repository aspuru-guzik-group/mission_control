import json
import os
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem

from . import e2e_flow_test_utils

from mc.a2g2.flow_generators.confgen_flow_generator import ConfgenFlowGenerator
from mc.a2g2.job_engines import a2g2_job_engine


def generate_test_mols():
    smiles_list = ['Cc1ccccc1']
    mols = []
    for smiles in smiles_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMultipleConfs(mol, numConfs=3)
        mol.smiles = smiles
        mols.append(mol)
    return mols

@unittest.skipUnless(*e2e_flow_test_utils.get_skip_args())
class ConfgenFlow_E2E_TestCase(e2e_flow_test_utils.E2E_Flow_BaseTestCase):
    def teardown_docker_env(self, docker_env=None):
        pass

    def get_flow_generator_classes(self):
        return [ConfgenFlowGenerator]

    def get_job_engine_class_spec(self):
        #return 'mc.a2g2.e2e_tests.test_confgen_flow_e2e.MockJobEngine'
        return None

    def generate_job_submission_cfg(self):  
        cfg = super().generate_job_submission_cfg()
        cfg.setdefault('env_vars', {})['CONFGEN_EXE'] = 'fooga'
        return cfg

    def test_flow(self):
        self.mc_client.flush_mc_db()
        self.a2g2_client.flush_a2g2_db()
        self.generate_molecule_library()
        self.create_flows()
        self.assertTrue(len(self.mc_client.fetch_tickable_flows()) > 0)
        try:
            self.run_flows_to_completion(max_ticks=15)
        except Exception as exception:
            jobs = self.mc_client.fetch_jobs()
            for job in jobs:
                error = job.get('error', '')
                print("job['error']:", error)
            raise exception
        self.assertTrue(self.combo_runner.tick_counter > 0)
        self.assert_domain_db_has_expected_state()

    def generate_molecule_library(self):
        mols = generate_test_mols()
        for mol in mols:
            self.a2g2_client.create_chemthing({'props': {'smiles': mol.smiles}})

    def create_flows(self):
        for chemthing in self.a2g2_client.query(
            q={'collection': 'chemthings'}
        )[:1]:
            self.create_confgen_flow(chemthing=chemthing)

    def create_confgen_flow(self, chemthing=None):
        flow_spec = {
            'flow_type': ConfgenFlowGenerator.flow_type,
            'input': {
                'smiles': chemthing['props']['smiles'],
                'confgen_params': {},
            },
        }
        flow = {'spec': json.dumps(flow_spec)}
        self.mc_client.create_flow(flow=flow)

    def assert_domain_db_has_expected_state(self):
        mols = generate_test_mols()
        mol_count = len(mols)
        conformer_count = 0
        for mol in mols: conformer_count += len(mol.GetConformers())
        expected_counts = {'ChemThing': mol_count + conformer_count }
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)


class MockJobEngine(object):
    def execute_command(self, command=None, job=None, cfg=None, **kwargs):
        use_real_engine = True
        if command == 'run_job_submission':
            if job['job_spec']['job_type'] == 'a2g2.jobs.confgen':
                use_real_engine = False
                raise NotImplementedError("implement fake confgen")
                #self.generate_fake_confgen_output(output_dir=output_dir)
        if use_real_engine:
            real_engine = a2g2_job_engine.A2G2JobEngine()
            real_engine.execute_command(command=command, job=job, cfg=cfg,
                                        **kwargs)

    def generate_fake_confgen_output(self, output_dir=None):
        mols = generate_test_mols()
        for mol_counter, mol in enumerate(mols):
            conformer_file_path = os.path.join(
                output_dir, 'mol_%s_conformers.sdf' % mol_counter)
            with open(conformer_file_path, 'w') as conformer_file:
                writer = Chem.SDWriter(conformer_file)
                for conformer in mol.GetConformers():
                    writer.write(mol, confId=conformer.GetId())
                writer.close()
