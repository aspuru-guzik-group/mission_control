import argparse
import gzip
import json
import os
import unittest

from mc.a2g2.flow_generators.qchem_flow_generator import QChemFlowGenerator
from mc.a2g2.job_modules.qchem.parse import tests as qchem_parse_tests

from . import e2e_flow_test_utils


class Fixtures(object):
    def __init__(self):
        self.molecule = {
            'charge': 0,
            'multiplicity': 1,
            'atoms': [
                {'element': 'C', 'x': i, 'y': i, 'z': i}
                for i in range(3)
            ]
        }
        self.precursors = {'precursor1': True}

    def generate_qchem_output_file_content(self):
        example_confgen_output_path = os.path.join(
            os.path.dirname(qchem_parse_tests.__file__), 'example_qchem.out.gz')
        return gzip.open(example_confgen_output_path, 'rb').read()

@unittest.skipUnless(*e2e_flow_test_utils.get_skip_args())
class QChemFlow_E2E_TestCase(e2e_flow_test_utils.E2E_Flow_BaseTestCase):
    def setUp(self):
        super().setUp()
        self.fixtures = Fixtures()
    
    def teardown_docker_env(self, docker_env=None): pass

    def get_flow_generator_classes(self): return [QChemFlowGenerator]

    def generate_job_submission_cfg(self):  
        cfg = super().generate_job_submission_cfg()
        cfg['a2g2.jobs.qchem.qchem'] = {
            'env_vars': {
                # fake confgen by calling this module' 
                'QCHEM_EXE': 'python -m mc.%s fake_qchem' % (__name__)
            }
        }
        return cfg

    def test_flow(self):
        try:
            self.mc_client.flush_mc_db()
            self.a2g2_client.flush_a2g2_db()
            self.create_flows()
            self.assertTrue(len(self.mc_client.fetch_tickable_flows()) > 0)
            self.run_flows_to_completion(timeout=20)
            self.assert_domain_db_has_expected_state()
        except Exception as exception:
            raise

    def create_flows(self):
        flow_spec = {
            'flow_type': QChemFlowGenerator.flow_type,
            'input': {
                'qchem_params': {
                    'molecule': self.fixtures.molecule,
                },
                'precursors': self.fixtures.precursors,
            },
        }
        flow = {'spec': json.dumps(flow_spec)}
        self.mc_client.create_flow(flow=flow)

    def assert_domain_db_has_expected_state(self):
        try:
            chemthings = self.a2g2_client.query(q={'collection': 'chemthings'})
            self.assert_has_expected_calc_chemthings(chemthings=chemthings)
        except:
            self.dump_db_state()
            raise

    def assert_has_expected_calc_chemthings(self, chemthings=None):
        actual_chemthings = self.filter_chemthings_by_type(
            chemthings=chemthings, _type='a2g2:type:calc')
        expected_partial_chemthings = [
            {
                'types': {'a2g2:type:calc': True},
                'precursors': self.fixtures.precursors,
            }
        ]
        self.assert_plucked_lists_equal(
            left=actual_chemthings,
            right=expected_partial_chemthings,
            keys=['types', 'precursors'])

#### Handlers for mocking executables ####

def handle_fake_qchem_command(args=None):
    # Simulate args for normal qchem call ('qchem infile oufile')
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    parsed_args, extra_args = parser.parse_known_args(args=args)
    fixtures = Fixtures()
    open(parsed_args.output_file,'wb').write(
        fixtures.generate_qchem_output_file_content())

if __name__ == '__main__':
    # Map fake calls to handlers.
    # Will be called in fake execution context, NOT in unittest context.
    parser = argparse.ArgumentParser(description='mocked executables')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    if parsed_args.command == 'fake_qchem':
        handle_fake_qchem_command(args=command_args)
