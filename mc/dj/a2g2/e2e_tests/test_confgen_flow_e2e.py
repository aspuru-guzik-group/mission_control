import json

from django.test import TestCase


class ConfgenFlow_E2E_TestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.a2g2_client = self.generate_a2g2_client()
        self.runner = self.generate_runner()

    def generate_a2g2_client(): raise NotImplementedError
    def generate_runner(): raise NotImplementedError

    def test_flow(self):
        self.generate_molecule_library()
        self.create_flows()
        self.run_flows_to_completion()
        self.assert_domain_db_has_expected_state()

    def generate_molecule_library(self):
        initial_smiles = []
        for smiles in initial_smiles:
            self.a2g2_client.create_molecule(smiles=smiles)

    def create_flows(self):
        for molecule in self.a2g2_client.query(q={'type': 'molecule'}):
            self.create_confgen_flow(molecule=molecule)

    def create_confgen_flow(self, molecule=None):
        flow_spec = {
            'type': 'confgen',
            'input': {
                'smiles': molecule['smiles'],
            },
        }
        flow = {'spec': json.dumps(flow_spec)}
        self.flow_client.create_flow(flow=flow)

    def run_flows_to_completion(self):
        tickable_flows = self.flow_client.get_tickable_flows()
        max_ticks = 30
        while len(tickable_flows) > 0:
            self.runner.tick()
            if self.runner.tick_counter > max_ticks:
                raise Exception("Exceeded max_ticks of '%s'" % max_ticks)

    def assert_domain_db_has_expected_state(self):
        expected_counts = {'ChemObj': 10}
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)
