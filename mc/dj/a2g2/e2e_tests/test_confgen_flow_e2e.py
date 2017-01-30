import json
from unittest.mock import MagicMock, Mock

from django.conf.urls import url, include
from django.test import TestCase, override_settings

from mc_utils import test_utils
from flow_engines.flow import Flow
from a2g2.a2g2_client.a2g2_client import A2G2_Client
from a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from a2g2.a2g2_dj import urls as _a2g2_dj_urls
from jobs import urls as _jobs_urls
from missions import urls as _missions_urls
from a2g2.tasks.flow import FlowTask
from a2g2.flow_generators.run_and_load import RunAndLoadFlowGenerator
from job_runners.action_processor import ActionProcessor


BASE_PATH = 'test_api'
BASE_URL = '/' + BASE_PATH
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_a2g2_dj_urls.__name__)),
    url(r'^%s' % BASE_PATH, include(_jobs_urls.__name__)),
    url(r'^%s' % BASE_PATH, include(_missions_urls.__name__)),
]

@override_settings(ROOT_URLCONF=__name__)
class ConfgenFlow_E2E_TestCase(TestCase):
    def setUp(self):
        super().setUp()
        test_utils.patch_request_client(request_client=self.client,
                                        json_methods=['get', 'post', 'patch'])
        self.a2g2_client = self.generate_a2g2_client()
        self.storage_client = self.generate_storage_client()
        self.action_processor = self.generate_action_processor()
        self.execution_client = self.generate_execution_client()
        self.job_dir_factory = MagicMock()
        self.flow_and_job_runner = self.generate_flow_and_job_runner()
        self.flow_client = self.flow_and_job_runner.flow_client
        self.flow_and_job_runner.tick = Mock(
            side_effect=self.flow_and_job_runner.tick)

    def generate_a2g2_client(self):
        return A2G2_Client(base_url=BASE_URL, request_client=self.client)

    def generate_storage_client(self):
        storage_client = Mock()
        storage_client.put.return_value = 'put_return'
        storage_client.get.return_value = 'get_return'
        return storage_client

    def generate_action_processor(self):
        action_processor = ActionProcessor()
        action_processor.register_handler(key='storage:put',
                                          handler=self.storage_client.put)
        action_processor.register_handler(key='storage:get',
                                          handler=self.storage_client.get)
        return action_processor

    def generate_execution_client(self):
        execution_client = MagicMock()
        execution_client.get_execution_state = Mock(return_value={})
        return execution_client

    def generate_flow_and_job_runner(self):
        return OdysseyPushRunner(
            flow_generator_classes=[ConfgenFlowGenerator],
            request_client=self.client,
            job_server_url=BASE_URL,
            flow_server_url=BASE_URL,
            job_dir_factory=self.job_dir_factory,
            job_runner_kwargs={'execution_client': self.execution_client,
                               'action_processor': self.action_processor,}
        )

    def test_flow(self):
        import logging
        logging.getLogger().setLevel(logging.DEBUG)

        self.generate_molecule_library()
        self.create_flows()
        self.assertTrue(len(self.flow_client.fetch_tickable_flows()) > 0)
        self.run_flows_to_completion()
        self.assertTrue(self.flow_and_job_runner.tick_counter > 0)
        self.assert_domain_db_has_expected_state()

    def generate_molecule_library(self):
        initial_smiles = ['smiles_%s' % i for i in range(3)]
        for smiles in initial_smiles:
            self.a2g2_client.create_mol({'props': {'smiles': smiles}})

    def create_flows(self):
        for mol in self.a2g2_client.query(q={'collection': 'mols'})[:1]:
            self.create_confgen_flow(mol=mol)

    def create_confgen_flow(self, mol=None):
        flow_spec = {
            'flow_type': ConfgenFlowGenerator.flow_type,
            'input': {
                'smiles': mol['props']['smiles'],
                'confgen_params': {},
            },
        }
        flow = {'spec': json.dumps(flow_spec)}
        self.flow_client.create_flow(flow=flow)

    def run_flows_to_completion(self):
        complete_flows = self.fetch_and_key_flows(
            query_params={'status': 'COMPLETED'})
        all_flows = self.fetch_and_key_flows()
        incomplete_flows = [flow for flow_uuid, flow in all_flows.items()
                            if flow_uuid not in complete_flows]
        #max_ticks = 30
        max_ticks = 4
        while len(incomplete_flows) > 0:
            if self.flow_and_job_runner.tick_counter > max_ticks:
                raise Exception("Exceeded max_ticks of '%s'" % max_ticks)
            self.flow_and_job_runner.tick()

    def fetch_and_key_flows(self, query_params=None):
        keyed_flows = {
            flow['uuid']: flow
            for flow in self.flow_client.fetch_flows(query_params=query_params)
        }
        return keyed_flows

    def assert_domain_db_has_expected_state(self):
        expected_counts = {'Mol': 10}
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)


class ConfgenFlowGenerator(object):
    flow_type = 'confgen'

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.add_task(
            key='confgen', 
            as_root=True,
            task=FlowTask(),
            input={
                'flow_spec': {
                    'flow_type': RunAndLoadFlowGenerator.flow_type,
                    'run_spec': {
                        'job_type': 'confgen',
                        'confgen': {
                            'smiles': flow_spec['input']['smiles'],
                            'params': flow_spec['input']['confgen_params'],
                        },
                        'post_exec_actions': [
                            {
                                'action': 'storage:put',
                                'params': {'src': '{ctx.completed_dir}'},
                                'output_to_ctx_target': 'data.output.raw_dir'
                            }
                        ]
                    },
                    'load_spec': {
                        'job_type': 'confgen:load',
                        'pre_build_actions': [
                            {
                                'action': 'storage:get',
                                'params': {
                                    'src': '{ctx.data.input.raw_dir}',
                                    'dest': '{ctx.job_dir}/raw_dir',
                                },
                                'output_to_ctx_target': 'data.input.raw_dir'
                            }
                        ]
                    }
                }
            })
        return flow

    @classmethod
    def get_dependencies(cls):
        return {
            'task_classes': set([FlowTask]),
            'flow_generator_classes': set([RunAndLoadFlowGenerator])
        }
