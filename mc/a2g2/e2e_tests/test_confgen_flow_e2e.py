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

class Fixtures(object):
    def __init__(self):
        self.mols = self.generate_mols()
        self.conformers = self.generate_conformers()
        self.confgen_params = self.generate_confgen_params()
    
    def generate_mols(self):
        smiles_list = ['Cc1ccccc1']
        return [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    def generate_conformers(self, num_confs=1):
        mol = self.mols[0]
        AllChem.EmbedMultipleConfs(mol, numConfs=num_confs)
        return mol.GetConformers()

    def generate_confgen_params(self):
        return {'param_%s' % i: 'value_%s' % i for i in range(3)}

    def conformer_to_xyz(self, conformer=None):
        return self.atoms_to_xyz(
            atoms=self.conformer_to_atoms(conformer),
            comment=self.generate_conformer_comment(conformer)
        )

    def atoms_to_xyz(self, atoms=None, comment=None):
        xyz_lines = []
        xyz_lines.append("%s" % len(atoms))
        xyz_lines.append("%s" % (comment or ""))
        for atom in atoms:
            xyz_lines.append("{element} {x:.4f} {y:.4f} {z:.4f}".format(**atom))
        xyz = "\n".join(xyz_lines)
        return xyz

    def conformer_to_atoms(self, conformer=None):
        atoms = []
        for i, atom in enumerate(conformer.GetOwningMol().GetAtoms()):
            pos = conformer.GetAtomPosition(i)
            coords = {coord: float("%.4f" % getattr(pos, coord))
                      for coord in ['x', 'y', 'z']}
            atoms.append({'element': atom.GetSymbol(), **coords})
        return atoms

    def generate_conformer_comment(self, conformer=None):
        return json.dumps({'rdkit_conformer_id': conformer.GetId()})

@unittest.skipUnless(*e2e_flow_test_utils.get_skip_args())
class ConfgenFlow_E2E_TestCase(e2e_flow_test_utils.E2E_Flow_BaseTestCase):
    def setUp(self):
        super().setUp()
        self.fixtures = Fixtures()
    
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
            self.run_flows_to_completion(timeout=20)
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
        for mol in self.fixtures.mols:
            self.a2g2_client.create_chemthing({
                'types': {'a2g2:type:mol2d': True},
                'props': {'a2g2:prop:smiles': Chem.MolToSmiles(mol)}
            })

    def create_flows(self):
        for chemthing in self.a2g2_client.query(
            q={'collection': 'chemthings'}
        )[:1]:
            self.create_confgen_flow(chemthing=chemthing)

    def create_confgen_flow(self, chemthing=None):
        flow_spec = {
            'flow_type': ConfgenFlowGenerator.flow_type,
            'input': {
                'smiles': chemthing['props']['a2g2:prop:smiles'],
                'confgen_params': {},
                'precursors': {
                    chemthing['uuid']: True,
                }
            },
        }
        flow = {'spec': json.dumps(flow_spec)}
        self.mc_client.create_flow(flow=flow)

    def assert_domain_db_has_expected_state(self):
        chemthings = self.a2g2_client.query(q={'collection': 'chemthings'})
        self.assert_has_expected_mol_chemthings(chemthings=chemthings)
        self.assert_has_expected_calc_chemthings(chemthings=chemthings)
        self.assert_has_expected_conformer_chemthings(chemthings=chemthings)

    def assert_has_expected_mol_chemthings(self, chemthings):
        actual_chemthings = self.filter_chemthings_by_type(
            chemthings=chemthings, _type='a2g2:type:mol2d')
        expected_chemthings = []
        for mol in self.fixtures.mols:
            expected_chemthings.append({
                'props': {
                    'a2g2:prop:smiles': Chem.MolToSmiles(mol),
                },
                'types': {'a2g2:type:mol2d': True}
            })
        self.assertEqual(len(actual_chemthings), len(expected_chemthings))
        def _sort_key(chemthing): return chemthing['props']['a2g2:prop:smiles']
        self.assert_plucked_lists_equal(
            left=sorted(actual_chemthings, key=_sort_key),
            right=sorted(expected_chemthings, key=_sort_key),
            keys=['props', 'types']
        )

    def filter_chemthings_by_type(self, chemthings=None, _type=None):
        return [chemthing for chemthing in chemthings
                if _type in chemthing['types']]

    def assert_plucked_lists_equal(self, left=None, right=None, keys=None):
        self.assertEqual([self._pluck(item, keys) for item in left],
                         [self._pluck(item, keys) for item in right])

    def _pluck(self, _dict, keys=None):
        return {k: v for k, v in _dict.items() if k in keys}

    def assert_has_expected_calc_chemthings(self, chemthings=None):
        actual_chemthings = self.filter_chemthings_by_type(
            chemthings=chemthings, _type='a2g2:type:calc:confgen')
        expected_chemthings = []
        for mol_chemthing in self.filter_chemthings_by_type(
            chemthings=chemthings, _type='a2g2:type:mol2d'):
            expected_chemthings.append({
                'types': {'a2g2:type:calc:confgen': True},
                'props': {
                    'a2g2:prop:confgen:parameters': (
                        self.fixtures.confgen_params)
                },
                'precursors': {mol_chemthing['uuid']: True}
            })
        self.assertEqual(len(actual_chemthings), len(expected_chemthings))
        def _sort_key(chemthing):
            return ','.join(sorted(chemthing['precursors'].keys()))
        self.assert_plucked_lists_equal(
            left=sorted(actual_chemthings, key=_sort_key),
            right=sorted(expected_chemthings, key=_sort_key),
            keys=['props', 'types', 'precursors']
        )

    def assert_has_expected_conformer_chemthings(self, chemthings=None):
        actual_chemthings = self.filter_chemthings_by_type(
            chemthings=chemthings, _type='a2g2:type:mol3d')
        expected_chemthings = []
        for calc_chemthing in self.filter_chemthings_by_type(
            chemthings=chemthings, _type='a2g2:type:calc:confgen'):
            for conformer in self.fixtures.conformers:
                expected_chemthings.append({
                    'types': {'a2g2:type:mol3d': True},
                    'props': {
                        'a2g2:prop:atoms': self.fixtures.conformer_to_atoms(
                            conformer),
                        'a2g2:prop:confgen:comment': (
                            self.fixtures.generate_conformer_comment(
                                conformer=conformer))
                    },
                    'precursors': {calc_chemthing['uuid']: True}
                })
        self.assertEqual(len(actual_chemthings), len(expected_chemthings))
        def _sort_key(chemthing):
            key_elements = [
                ";".join(sorted(chemthing['precursors'].keys())),
                chemthing['props']['a2g2:prop:confgen:comment']
            ]
            return ','.join(key_elements)
        self.assert_plucked_lists_equal(
            left=sorted(actual_chemthings, key=_sort_key),
            right=sorted(expected_chemthings, key=_sort_key),
            keys=['props', 'types', 'precursors']
        )

#### Handlers for mocking executables ####

def handle_confgen_command(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir')
    parsed_args, extra_args = parser.parse_known_args(args=args)
    outdir = parsed_args.outdir
    xyz_dir = os.path.join(outdir, 'conformers')
    os.makedirs(xyz_dir, exist_ok=True)
    fixtures = Fixtures()
    for i, conformer in enumerate(fixtures.conformers):
        xyz_path = os.path.join(xyz_dir, 'conformer_%s.xyz' % i)
        open(xyz_path, 'w').write(fixtures.conformer_to_xyz(conformer))
    from mc.a2g2.job_modules.confgen.constants import CONFGEN_PARAMS_FILENAME
    confgen_params_outfile_path = os.path.join(outdir, CONFGEN_PARAMS_FILENAME)
    with open(confgen_params_outfile_path, 'w') as f:
        f.write(json.dumps(fixtures.confgen_params))

if __name__ == '__main__':
    # Logic to fake calls to executables.
    # This will be called in the fake execution context, NOT in the
    # unittest context.
    parser = argparse.ArgumentParser(description='mocked executables')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    if parsed_args.command == 'fake_confgen':
        handle_confgen_command(args=command_args)
