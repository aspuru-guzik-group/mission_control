import argparse
import json
import unittest

from . import e2e_flow_test_utils
from .fixtures.confgen import ConfgenFixtures

from mc.a2g2.flow_generators.confgen_flow_generator import ConfgenFlowGenerator


@unittest.skipUnless(*e2e_flow_test_utils.get_skip_args())
class ConfgenFlow_E2E_TestCase(e2e_flow_test_utils.E2E_Flow_BaseTestCase):
    def setUp(self):
        super().setUp()
        self.confgen_fixtures = ConfgenFixtures()
    
    def teardown_docker_env(self, docker_env=None): pass

    def get_flow_generator_classes(self): return [ConfgenFlowGenerator]

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

    def generate_molecule_library(self):
        for mol in self.confgen_fixtures.mols:
            self.a2g2_client.create_chemthing({
                'types': {'a2g2:type:mol2d': True},
                'props': {'a2g2:prop:smiles': (
                    self.confgen_fixtures.Chem.MolToSmiles(mol))}
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
        for mol in self.confgen_fixtures.mols:
            expected_chemthings.append({
                'props': {
                    'a2g2:prop:smiles': self.confgen_fixtures.Chem.MolToSmiles(
                        mol),
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
                        self.confgen_fixtures.confgen_params)
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
            for conformer in self.confgen_fixtures.conformers:
                expected_chemthings.append({
                    'types': {'a2g2:type:mol3d': True},
                    'props': {
                        'a2g2:prop:atoms': (
                            self.confgen_fixtures.conformer_to_atoms(conformer)
                        ),
                        'a2g2:prop:confgen:comment': (
                            self.confgen_fixtures.generate_conformer_comment(
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
if __name__ == '__main__':
    # Map fake calls to handlers.
    # Will be called in fake execution context, NOT in unittest context.
    parser = argparse.ArgumentParser(description='mocked executables')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    if parsed_args.command == 'fake_confgen':
        ConfgenFixtures.handle_confgen_command(args=command_args)
