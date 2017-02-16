import json
import sys
import unittest

import pexpect

from flow_engines.flow import Flow
from a2g2.a2g2_client.a2g2_client import A2G2_Client
from a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from storage_client.storage_client import MissionControlStorageClient
from a2g2.task_engines.job_task_engine import JobTaskEngine
from job_runners.action_processor import ActionProcessor

from .docker_utils import DockerEnv
from . import storage_action_handlers
from job_runners.ssh_control_socket_client import SSHControlSocketClient


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
                                    'data.output.storage_meta')
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
                                    '.storage_meta}}'
                                ),
                            },
                            'target': 'task.input.job_spec.input.storage_meta',
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
                                    'storage_meta': {
                                        'template': (
                                            '{{ctx.job_spec.input.storage_meta}}'
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
        self.job_submission_factory = self.generate_job_submission_factory()
        self.mc_runner = self.generate_mc_runner()
        self.mc_client = self.mc_runner.mc_client

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

    def generate_job_submission_factory(self):
        a2g2_client_cfg_json = json.dumps({
            'base_url': self.a2g2_client.base_url
        })

        class JobSubmissionFactory(object):
            def build_job_submission(self, job=None):
                job_type = job['job_spec']['job_type']
                if job_type == 'confgen':
                    from a2g2.job_dir_builders.confgen.confgen \
                            import ConfgenJobDirBuilder
                    job_dir_meta = ConfgenJobDirBuilder.build_odyssey_dir(
                        job=job
                    )
                elif job_type == 'confgen:load':
                    from a2g2.job_dir_builders.confgen_load.confgen_load \
                            import ConfgenLoadJobDirBuilder
                    job_dir_meta = ConfgenLoadJobDirBuilder.build_odyssey_dir(
                        job=job,
                        a2g2_client_cfg_json=a2g2_client_cfg_json
                    )
                else:
                    raise Exception("Unknown job type '%s',"
                                    " can't build job dir" % job_type)
                submission_meta = job_dir_meta
                return submission_meta
        return JobSubmissionFactory()

    def generate_mc_runner(self):
        return OdysseyPushRunner(
            ssh_client=self.generate_connected_ssh_client(),
            action_processor=self.action_processor,
            flow_generator_classes=[ConfgenFlowGenerator],
            mc_server_url='http://' + self.docker_env.mc_ip_addr + '/missions/',
            job_submission_factory=self.job_submission_factory,
            job_runner_kwargs={'action_processor': self.action_processor}
        )

    def generate_connected_ssh_client(self):
        def auth_fn(*args, ssh_client=None, ssh_cmd=None, **kwargs):
            child = pexpect.spawn(ssh_cmd, timeout=5)
            child.expect('.+')
            child.sendline(self.docker_env.odyssey_user_password)
            child.expect('.+')
            child.sendline(ssh_client.user)
            child.wait()
            if child.exitstatus != 0:
                raise Exception("Authorization failed")
            print("AUTHORIZED", file=sys.stderr)

        ssh_client = SSHControlSocketClient(
            user=self.docker_env.odyssey_user_username,
            host=self.docker_env.odyssey_ip_addr,
            auth_fn=auth_fn,
        )
        ssh_client.connect()
        return ssh_client

    def test_flow(self):
        self.mc_client.flush_mc_db()
        self.a2g2_client.flush_a2g2_db()
        self.generate_molecule_library()
        self.create_flows()
        self.assertTrue(len(self.mc_client.fetch_tickable_flows()) > 0)
        try:
            self.run_flows_to_completion()
            err = None
        except Exception as error:
            err = error
        jobs = self.mc_client.fetch_jobs()
        print("len: ", len(jobs))
        print(json.dumps(jobs, indent=2))
        if err: raise err
        self.assertTrue(self.mc_runner.tick_counter > 0)
        self.assert_domain_db_has_expected_state()

    def generate_molecule_library(self):
        initial_smiles = ['smiles_%s' % i for i in range(1)]
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
        self.mc_client.create_flow(flow=flow)

    def run_flows_to_completion(self):
        max_ticks = 20
        incomplete_flows = self.get_incomplete_flows()
        while len(incomplete_flows) > 0:
            self.mc_runner.tick()
            if self.mc_runner.tick_counter > max_ticks:
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
            for flow in self.mc_client.fetch_flows(query_params=query_params)
        }
        return keyed_flows

    def assert_domain_db_has_expected_state(self):
        expected_counts = {'ChemThing': 10}
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)

