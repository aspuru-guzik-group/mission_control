import json
import os
import time
import unittest

import pexpect
from rdkit import Chem
from rdkit.Chem import AllChem

from mc.flow_engines.flow import Flow
from mc.a2g2.a2g2_client.a2g2_client import A2G2_Client
from mc.a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from mc.storage_client.storage_client import MissionControlStorageClient
from mc.a2g2.task_engines.job_task_engine import JobTaskEngine
from mc.job_runners.action_processor import ActionProcessor
from mc.job_runners.ssh_control_socket_client import SSHControlSocketClient
from mc.a2g2.job_dir_builders.a2g2_job_engine.a2g2_job_engine_dir_builder import (
    A2G2JobEngineDirBuilder)
from mc.a2g2.job_engines import a2g2_job_engine
from mc.a2g2.action_handlers import storage_action_handlers
from .docker.docker_utils import DockerEnv


def generate_test_mols():
    smiles_list = ['Cc1ccccc1']
    mols = []
    for smiles in smiles_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMultipleConfs(mol, numConfs=3)
        mol.smiles = smiles
        mols.append(mol)
    return mols

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
                        'job_type': 'confgen_load',
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
                                        'template': (
                                            '{{ctx.job_dir}}/dir_to_parse')
                                    }
                                },
                            },
                            {
                                'action': 'set_ctx_value',
                                'description': 'set name of dir_to_parse',
                                'params': {
                                    'value': 'dir_to_parse',
                                    'target': 'data.input.dir_to_parse',
                                }
                            },
                        ]
                    }
                },
            }
        )
        return flow

    @classmethod
    def get_dependencies(cls):
        return {
            'task_engines': set([JobTaskEngine()]),
        }

@unittest.skipUnless(os.environ.get('MC_RUN_DOCKER_TESTS', False),
                     'skipping because MC_RUN_DOCKER_TESTS is not set')
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
            def build_job_submission(self, job=None, output_dir=None):
                job_dir_meta = A2G2JobEngineDirBuilder.build_odyssey_dir(
                    job=job,
                    cfg={
                        'A2G2_CLIENT_CFG_JSON': a2g2_client_cfg_json,
                        'A2G2_JOB_ENGINE_CLASS': (
                            'mc.a2g2.e2e_tests.'
                            'test_confgen_flow_e2e.MockJobEngine'
                        )
                    },
                    output_dir=output_dir
                )
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
            if child.exitstatus != 0: raise Exception("Authorization failed")

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
        except Exception as error:
            jobs = self.mc_client.fetch_jobs()
            print(json.dumps(jobs, indent=2))
            raise error
        self.assertTrue(self.mc_runner.tick_counter > 0)
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

    def run_flows_to_completion(self):
        max_ticks = 10
        tick_interval = 1
        incomplete_flows = self.get_incomplete_flows()
        while len(incomplete_flows) > 0:
            self.mc_runner.tick()
            if self.mc_runner.tick_counter > max_ticks:
                raise Exception("Exceeded max_ticks of '%s'" % max_ticks)
            incomplete_flows = self.get_incomplete_flows()
            time.sleep(tick_interval)

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
        mols = generate_test_mols()
        mol_count = len(mols)
        conformer_count = 0
        for mol in mols: conformer_count += len(mol.GetConformers())
        expected_counts = {'ChemThing': mol_count + conformer_count }
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)


class MockJobEngine(object):
    def execute_job(self, job=None, cfg=None, output_dir=None, ctx_dir=None):
        job_type = job['job_spec']['job_type']
        if job_type == 'confgen':
            self.generate_fake_confgen_output(output_dir=output_dir)
        else:
            real_engine = a2g2_job_engine.A2G2JobEngine()
            real_engine.execute_job(job=job, cfg=cfg, output_dir=output_dir,
                                    ctx_dir=ctx_dir)

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
