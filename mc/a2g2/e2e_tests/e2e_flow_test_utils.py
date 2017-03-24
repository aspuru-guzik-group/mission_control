import json
import os
import time
import unittest

import pexpect

from mc.a2g2.a2g2_client.a2g2_client import A2G2_Client
from mc.a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from mc.storage_client.storage_client import MissionControlStorageClient
from mc.execution_clients.ssh_control_socket_client import (
    SSHControlSocketClient)
from mc.a2g2.job_dir_builders.a2g2_job_engine.a2g2_job_engine_dir_builder import (
    A2G2JobEngineDirBuilder)
from .docker.docker_utils import DockerEnv

def get_skip_args():
    return (os.environ.get('MC_RUN_DOCKER_TESTS', False),
            'skipping because MC_RUN_DOCKER_TESTS is not set')

class E2E_Flow_BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.docker_env = self.setup_docker_env()
        self.a2g2_client = self.generate_a2g2_client()
        self.storage_client = self.generate_storage_client()
        self.job_submission_factory = self.generate_job_submission_factory()
        self.task_handler = self.generate_task_handler()
        self.mc_runner = self.generate_mc_runner()
        self.mc_client = self.mc_runner.mc_client

    def tearDown(self):
        self.teardown_docker_env(docker_env=self.docker_env)

    def setup_docker_env(self):
        docker_env = DockerEnv()
        docker_env.setup()
        return docker_env

    def teardown_docker_env(self, docker_env=None):
        docker_env.teardown()

    def generate_a2g2_client(self):
        return A2G2_Client(
            base_url='http://' + self.docker_env.mc_ip_addr + '/a2g2/')

    def generate_storage_client(self):
        storage_client = MissionControlStorageClient(
            base_url='http://' + self.docker_env.mc_ip_addr + '/storage/')
        return storage_client

    def generate_task_handler(self):
        class TaskHandler(object):
            def __init__(self):
                self.handlers = {}

            def tick_task(self, *args, task=None, task_context=None, **kwargs):
                handler = self.handlers[task['task_type']]
                handler.tick_task(*args, task=None, task_context=None, **kwargs)

        return TaskHandler()

    def _upload_wrapper(self, *args, params=None, ctx=None, **kwargs):
        pass

    def _download_wrapper(self, *args, params=None, ctx=None, **kwargs):
        pass

    def generate_job_submission_factory(self):
        a2g2_client_cfg_json = json.dumps({
            'base_url': self.a2g2_client.base_url
        })

        job_engine_class_spec = self.get_job_engine_class_spec()

        class JobSubmissionFactory(object):
            def build_job_submission(self, job=None, output_dir=None):
                job_dir_meta = A2G2JobEngineDirBuilder.build_odyssey_dir(
                    job=job,
                    cfg={
                        'A2G2_CLIENT_CFG_JSON': a2g2_client_cfg_json,
                        'A2G2_JOB_ENGINE_CLASS': job_engine_class_spec,
                    },
                    output_dir=output_dir
                )
                submission_meta = job_dir_meta
                return submission_meta
        return JobSubmissionFactory()

    def get_job_engine_class_spec(self): raise NotImplementedError

    def generate_mc_runner(self):
        return OdysseyPushRunner(
            ssh_client=self.generate_connected_ssh_client(),
            task_handler=self.task_handler,
            flow_generator_classes=self.get_flow_generator_classes(),
            mc_server_url='http://' + self.docker_env.mc_ip_addr + '/missions/',
            job_submission_factory=self.job_submission_factory,
        )

    def get_flow_generator_classes(self): raise NotImplementedError

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

    def run_flows_to_completion(self, max_ticks=10, tick_interval=1):
        incomplete_flows = self.get_incomplete_flows()
        while len(incomplete_flows) > 0:
            self.mc_runner.tick()
            if self.mc_runner.tick_counter > max_ticks:
                error_msg = (
                    "Exceeded max_ticks of '{max_ticks}.'\njobs:\n{jobs}"
                ).format(
                    max_ticks=max_ticks,
                    jobs=json.dumps(self.mc_client.fetch_jobs(), indent=2)
                )
                raise Exception(error_msg)
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
