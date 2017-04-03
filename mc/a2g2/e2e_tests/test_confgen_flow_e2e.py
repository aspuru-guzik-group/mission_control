import argparse
import json
import os
import textwrap
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
import yaml

from . import e2e_flow_test_utils

from mc.a2g2.flow_generators.confgen_flow_generator import ConfgenFlowGenerator


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
        cfg['a2g2.jobs.confgen.confgen'] = {
            'env_vars': {
                # fake confgen by calling this module' 
                'CONFGEN_EXE': 'python -m mc.%s fake_confgen' % (__name__)
            }
        }
        return cfg

    def test_flow(self):
        try:
            self.mc_client.flush_mc_db()
            self.a2g2_client.flush_a2g2_db()
            self.generate_molecule_library()
            self.create_flows()
            self.assertTrue(len(self.mc_client.fetch_tickable_flows()) > 0)
            self.run_flows_to_completion(timeout=15)
            self.assertTrue(self.combo_runner.tick_counter > 0)
            self.assert_domain_db_has_expected_state()
        except Exception as exception:
            #self.dump_db_state()
            raise exception

    def dump_db_state(self):
        jobs = self.mc_client.fetch_jobs()
        for i, job in enumerate(jobs):
            print("job %s:\n" % i, self.dump_obj(job))
        flows = self.mc_client.fetch_flows()
        for i, flow in enumerate(flows):
            print("flow %s:\n" % i, self.dump_obj(flow))

    def dump_obj(self, obj):
        return textwrap.indent(yaml.dump(obj), prefix=' ' * 2)

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

#### Handlers for mocking executables ####

def handle_confgen_command(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir')
    parsed_args, extra_args = parser.parse_known_args(args=args)
    outdir = parsed_args.outdir
    mols = generate_test_mols()

    def _conformer_to_xyz_str(conformer):
        return _atoms_to_xyz_str(_conformer_to_atoms(conformer))

    def _conformer_to_atoms(conformer):
        atoms = []
        for i, atom in enumerate(conformer.GetOwningMol().GetAtoms()):
            pos = conformer.GetAtomPosition(i)
            atoms.append([atom.GetSymbol(), [pos.x, pos.y, pos.z]])
        return atoms

    def _atoms_to_xyz_str(atoms=None, comment=None):
        xyz_lines = []
        xyz_lines.append("%s" % len(atoms))
        xyz_lines.append("%s" % (comment or ""))
        for atom in atoms:
            xyz_lines.append("%s %.4f %.4f %.4f" % (atom[0], *atom[1]))
        xyz_str = "\n".join(xyz_lines)
        return xyz_str

    xyz_dir = os.path.join(outdir, 'conformers')
    os.makedirs(xyz_dir, exist_ok=True)
    for mol_counter, mol in enumerate(mols):
        for i, conformer in enumerate(mol.GetConformers()):
            xyz_path = os.path.join(xyz_dir, 'conformer_%s.xyz' % i)
            open(xyz_path, 'w').write(_conformer_to_xyz_str(conformer))

if __name__ == '__main__':
    # Logic to fake calls to executables.
    # This will be called in the fake execution context, NOT in the
    # unittest context.
    parser = argparse.ArgumentParser(description='mocked executables')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    if parsed_args.command == 'fake_confgen':
        handle_confgen_command(args=command_args)
