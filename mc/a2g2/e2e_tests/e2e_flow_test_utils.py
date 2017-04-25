import json
import logging
import operator
import os
import sys
import tempfile
import textwrap
import time
import unittest

import yaml
import pexpect

from mc.a2g2.a2g2_client.a2g2_client import A2G2_Client
from mc.a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from mc.storage_client.storage_client import MissionControlStorageClient
from mc.execution_clients.ssh_control_socket_client import (
    SSHControlSocketClient)
from mc.a2g2.job_engines import a2g2_job_engine
from mc.a2g2.task_handlers.a2g2_task_handler import A2G2TaskHandler

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
        self.task_handler = self.generate_task_handler()
        self.combo_runner = self.generate_combo_runner()
        self.mc_client = self.combo_runner.mc_client

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
            def tick_task(self, *args, task=None, **kwargs):
                try:
                    A2G2TaskHandler.tick_task(task=task, **kwargs)
                except Exception as exception:
                    msg = ("Could not handle task with type"
                           " '%s'" % task['task_type'])
                    logging.exception(msg)
                    raise

            def compile_tasks(self, *args, **kwargs):
                return A2G2TaskHandler.compile_tasks(*args, **kwargs)

        return TaskHandler()

    def generate_combo_runner(self):
        return OdysseyPushRunner(
            task_handler=self.task_handler,
            job_submission_factory=self.generate_job_submission_factory(),
            ssh_client=self.generate_connected_ssh_client(),
            flow_generator_classes=self.get_flow_generator_classes(),
            mc_server_url='http://' + self.docker_env.mc_ip_addr + '/missions/',
        )

    def generate_job_submission_factory(self):
        cfg = self.generate_job_submission_cfg()
        class JobSubmissionFactory(object):
            def build_job_submission(self, job=None, submission_dir=None):
                submission_dir = submission_dir or tempfile.mkdtemp(
                    prefix='sf.')
                job_engine = a2g2_job_engine.A2G2JobEngine()
                submission_meta = job_engine.build_job_submission(
                    job=job,
                    cfg=cfg,
                    submission_dir=submission_dir
                )
                return submission_meta
        return JobSubmissionFactory()

    def generate_job_submission_cfg(self):  
        cfg = {
            'a2g2_client_cfg_json': json.dumps({
                'base_url': self.a2g2_client.base_url
            }),
            'job_engine': {
                'engine_class_spec': self.get_job_engine_class_spec(),
                'entrypoint_preamble': self.get_entrypoint_preamble(),
                'job_engine_exe': 'python -m %s' % (a2g2_job_engine.__name__),
            }
        }
        return cfg

    def get_entrypoint_preamble(self):
        cluster_conda_env_root = self.docker_env.compose_cfg['services']\
                ['odyssey']['environment']['ODYSSEY_A2G2_CONDA_ENV_PATH']
        preamble = textwrap.dedent(
            '''
            source {conda_env_root}/bin/activate {conda_env_root}
            '''
        ).strip().format(conda_env_root=cluster_conda_env_root)
        return preamble

    def get_job_engine_class_spec(self): return None

    def get_flow_generator_classes(self): return []

    def generate_connected_ssh_client(self):
        def auth_fn(*args, ssh_client=None, ssh_cmd=None, **kwargs):
            child = pexpect.spawn(ssh_cmd, timeout=5)
            child.expect('Password:')
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

    def run_flows_to_completion(self, max_ticks=50, tick_interval=0.5,
                                timeout=30, job_queue=None):
        start_time = time.time()
        self.combo_runner.job_queue_key = job_queue['uuid']
        incomplete_flows = self.get_incomplete_flows()
        while len(incomplete_flows) > 0:
            self.combo_runner.tick()
            error_msg = None
            if self.combo_runner.tick_counter > max_ticks:
                error_msg = ("Exceeded max_ticks of '{max_ticks}.'").format(
                    max_ticks=max_ticks)
            elif (time.time() - start_time) > timeout:
                error_msg = ("Exceeded timeout of {timeout} seconds.").format(
                    timeout=timeout)
            if error_msg: raise Exception(error_msg)
            incomplete_flows = self.get_incomplete_flows()
            time.sleep(tick_interval)

    def get_incomplete_flows(self):
        all_flows = self.fetch_and_key_flows()
        incomplete_flows = [flow for flow_uuid, flow in all_flows.items()
                            if flow['status'] not in ['COMPLETED', 'FAILED']]
        return incomplete_flows

    def fetch_and_key_flows(self, query_params=None):
        keyed_flows = {
            flow['uuid']: flow
            for flow in self.mc_client.fetch_flows(query_params=query_params)
        }
        return keyed_flows

    def dump_db_state(self, stream=sys.stdout):
        self.dump_db_jobs(stream=stream)
        self.dump_db_flows(stream=stream)

    def dump_db_jobs(self, stream=sys.stdout):
        jobs = self.mc_client.fetch_jobs()
        self.dump_db_objects(db_objects=jobs, label='job', stream=stream)

    def dump_db_objects(self, db_objects=None, label='', stream=sys.stdout):
        sorted_objects = sorted(db_objects, key=operator.itemgetter('created'))
        for i, obj in enumerate(sorted_objects):
            stream.write("{label} {i}:\n".format(label=label, i=i))
            stream.write(self.obj_to_yaml(obj))

    def obj_to_yaml(self, obj):
        return textwrap.indent(yaml.dump(obj), prefix=' ' * 2)

    def dump_db_flows(self, stream=sys.stdout):
        flows = self.mc_client.fetch_flows()
        self.dump_db_objects(db_objects=flows, label='flow', stream=stream)

    def filter_chemthings_by_type(self, chemthings=None, _type=None):
        return [chemthing for chemthing in chemthings
                if _type in chemthing['types']]

    def assert_plucked_lists_equal(self, left=None, right=None, keys=None):
        self.assertEqual([self._pluck(item, keys) for item in left],
                         [self._pluck(item, keys) for item in right])

    def _pluck(self, _dict, keys=None):
        return {k: v for k, v in _dict.items() if k in keys}
