import json
from unittest.mock import Mock

from djanog.conf.urls import url, include
from django.test import TestCase, override_settings

from mc import urls as _mc_urls
from mc_utils import test_utils
from a2g2.a2g2_client import A2G2_Client
from a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)

urlpatterns = [
    url(r'^', include(_mc_urls.__name__))
]

@override_settings(ROOT_URLCONF=__name__)
class ConfgenFlow_E2E_TestCase(TestCase):
    def setUp(self):
        super().setUp()
        test_utils.patch_request_client(request_client=self.client,
                                        json_methods=['get', 'post', 'patch'])
        self.a2g2_client = self.generate_a2g2_client()
        self.execution_client = Mock()
        self.flow_and_job_runner = self.generate_flow_and_job_runner()

    def generate_a2g2_client(self):
        self.a2g2_client = A2G2_Client(base_url='/', request_client=self.client)

    def generate_flow_and_job_runner(self):
        self.flow_and_job_runner = OdysseyPushRunner(
            request_client=self.request_client,
            job_server_url='/',
            flow_server_url='/',
            job_dir_factory=self.job_dir_factory,
            job_runner_kwargs={'execution_client': self.execution_client,
                               'transfer_client': self.transfer_client},
        )

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
