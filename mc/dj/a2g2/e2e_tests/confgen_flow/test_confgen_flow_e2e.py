import json
import tempfile
import unittest
from unittest.mock import MagicMock, Mock

from flow_engines.flow import Flow
from a2g2.a2g2_client.a2g2_client import A2G2_Client
from a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from storage_client.storage_client import MissionControlStorageClient
from a2g2.task_engines.job_task_engine import JobTaskEngine
from job_runners.action_processor import ActionProcessor

from .docker_utils import DockerEnv
from . import storage_action_handlers


class ConfgenFlowGenerator(object):
    flow_type = 'confgen'

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.data['flow_spec'] = flow_spec
        flow.add_task(
            key='confgen_run', 
            as_root=True,
            task={
                'task_engine': JobTaskEngine.__name__,
                'input': {
                    'job_spec': {
                        'job_type': 'confgen',
                        'confgen': {
                            'smiles': flow_spec['input']['smiles'],
                            'params': flow_spec['input']['confgen_params'],
                        },
                        'post_exec_actions': [
                            {
                                'action': 'storage:upload',
                                'params': {
                                    'src': {
                                        'template': '{{ctx.completed_dir}}'
                                    }
                                },
                                'output_to_ctx_target': (
                                    'data.output.raw_dir_storage_params')
                            }
                        ]
                    }
                },
                'status': 'PENDING',
            }
        )
        flow.add_task(
            key='confgen_load',
            precursor_keys=['confgen_run'],
            task={
                'pre_start_actions': [
                    {
                        'action': 'set_ctx_value',
                        'description': 'wire output from job task to job input',
                        'params': {
                            'value': {
                                'template': (
                                    '{{ctx.flow.tasks.confgen_run.output'
                                    '.raw_dir_storage_params}}'
                                ),
                            },
                            'target': ('task.input.job_spec.input'
                                       '.json_storage_params'),
                        }
                    }
                ],
                'task_engine': JobTaskEngine.__name__,
                'input': {
                    'job_spec': {
                        'job_type': 'confgen:load',
                        'pre_build_actions': [
                            {
                                'action': 'storage:download',
                                'params': {
                                    'json_src_params': {
                                        'template': (
                                            '{{ctx.job_spec.input.'
                                            'json_storage_params}}'
                                        ),
                                    },
                                    'dest': {
                                        'template': '{{ctx.job_dir}}/raw_dir'
                                    }
                                },
                                'output_to_ctx_target': 'data.input.raw_dir'
                            }
                        ]
                    }
                },
                'status': 'PENDING'
            }
        )
        return flow

    @classmethod
    def get_dependencies(cls):
        return {
            'task_engines': set([JobTaskEngine()]),
        }

class ConfgenFlow_E2E_TestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.docker_env = self.setup_docker_env()
        self.a2g2_client = self.generate_a2g2_client()
        self.storage_client = self.generate_storage_client()
        self.action_processor = self.generate_action_processor()
        self.execution_client = self.generate_execution_client()
        self.job_dir_factory = MagicMock()
        self.flow_and_job_runner = self.generate_flow_and_job_runner()
        self.flow_client = self.flow_and_job_runner.flow_client
        self.flow_and_job_runner.tick = Mock(
            side_effect=self.flow_and_job_runner.tick)

    def tearDown(self):
        #self.docker_env.teardown()
        pass

    def setup_docker_env(self):
        docker_env = DockerEnv()
        docker_env.setup()
        return docker_env

    def generate_a2g2_client(self):
        return A2G2_Client(
            base_url='http://' + self.docker_env.mc_ip_addr + '/a2g2/')

    def generate_storage_client(self):
        storage_client = MissionControlStorageClient(
            base_url='http://' + self.docker_env.mc_ip_addr + '/storage/')
        return storage_client

    def generate_action_processor(self):
        action_processor = ActionProcessor()
        action_processor.register_handler(key='storage:upload',
                                          handler=self._upload_wrapper)
        action_processor.register_handler(key='storage:download',
                                          handler=self._download_wrapper)
        return action_processor

    def _upload_wrapper(self, *args, params=None, ctx=None, **kwargs):
        return storage_action_handlers.upload_action_handler(
            self.storage_client, *args, params=params, ctx=ctx)

    def _download_wrapper(self, *args, params=None, ctx=None, **kwargs):
        return storage_action_handlers.download_action_handler(
            self.storage_client, *args, params=params, ctx=ctx)

    def generate_execution_client(self):
        execution_client = MagicMock()
        execution_client.get_execution_state = Mock(return_value={})
        def mock_start_execution(job=None):
            output_dir = tempfile.mkdtemp()
            if job['job_spec']['job_type'] == 'confgen':
                self.populate_confgen_dir(output_dir=output_dir)
            if job['job_spec']['job_type'] == 'confgen:load': pass
            return {'dir': output_dir}
        execution_client.start_execution.side_effect = mock_start_execution
        return execution_client

    def populate_confgen_dir(self, output_dir=None):
        pass

    def generate_flow_and_job_runner(self):
        return OdysseyPushRunner(
            action_processor=self.action_processor,
            flow_generator_classes=[ConfgenFlowGenerator],
            job_server_url='http://' + self.docker_env.mc_ip_addr + '/jobs/',
            flow_server_url=('http://' + self.docker_env.mc_ip_addr 
                             + '/missions/'),
            job_dir_factory=self.job_dir_factory,
            job_runner_kwargs={'execution_client': self.execution_client,
                               'action_processor': self.action_processor,}
        )

    def test_flow(self):
        self.generate_molecule_library()
        self.create_flows()
        self.assertTrue(len(self.flow_client.fetch_tickable_flows()) > 0)
        self.run_flows_to_completion()
        self.assertTrue(self.flow_and_job_runner.tick_counter > 0)
        self.assert_domain_db_has_expected_state()

    def generate_molecule_library(self):
        initial_smiles = ['smiles_%s' % i for i in range(3)]
        for smiles in initial_smiles:
            self.a2g2_client.create_chemthing({'props': {'smiles': smiles}})

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
        self.flow_client.create_flow(flow=flow)

    def run_flows_to_completion(self):
        max_ticks = 20
        incomplete_flows = self.get_incomplete_flows()
        while len(incomplete_flows) > 0:
            self.flow_and_job_runner.tick()
            if self.flow_and_job_runner.tick_counter > max_ticks:
                raise Exception("Exceeded max_ticks of '%s'" % max_ticks)
            incomplete_flows = self.get_incomplete_flows()

    def get_incomplete_flows(self):
        complete_flows = self.fetch_and_key_flows(
            query_params={'status': 'COMPLETED'})
        all_flows = self.fetch_and_key_flows()
        incomplete_flows = [flow for flow_uuid, flow in all_flows.items()
                            if flow_uuid not in complete_flows]
        return incomplete_flows

    def fetch_and_key_flows(self, query_params=None):
        keyed_flows = {
            flow['uuid']: flow
            for flow in self.flow_client.fetch_flows(query_params=query_params)
        }
        return keyed_flows

    def assert_domain_db_has_expected_state(self):
        expected_counts = {'ChemThing': 10}
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)

