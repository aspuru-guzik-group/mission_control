import argparse
import unittest

from mc.flow_engines.flow_engine import FlowEngine
from mc.a2g2.flow_generators import reaxys_flow_spec_generator

from . import e2e_flow_test_utils
from .fixtures.confgen import ConfgenFixtures
from .fixtures.qchem import QChemFixtures
from .fixtures.a2g2_db import A2G2_DB_Fixtures

import logging
logging.basicConfig(level=logging.INFO)


@unittest.skipUnless(*e2e_flow_test_utils.get_skip_args())
class QChemFlow_E2E_TestCase(e2e_flow_test_utils.E2E_Flow_BaseTestCase):
    def setUp(self):
        super().setUp()
        self.qchem_fixtures = QChemFixtures()
    
    def teardown_docker_env(self, docker_env=None): pass

    def generate_job_submission_cfg(self):  
        cfg = super().generate_job_submission_cfg()
        cfg['a2g2.jobs.confgen.confgen'] = {
            'env_vars': {
                # fake confgen by calling this module' 
                'CONFGEN_EXE': 'python -m mc.%s fake_confgen' % (__name__)
            }
        }
        cfg['a2g2.jobs.qchem.qchem'] = {
            'env_vars': {
                # fake confgen by calling this module' 
                'QCHEM_EXE': 'python -m mc.%s fake_qchem' % (__name__)
            }
        }
        cfg['a2g2.jobs.a2g2_db.query'] = {
            'env_vars': {
                # fake confgen by calling this module' 
                'A2G2_DB_EXE': 'python -m mc.%s fake_a2g2_db' % (__name__)
            }
        }
        return cfg

    def get_job_engine_class_spec(self):
        return 'mc.%s.FakeJobEngine' % (__name__)

    def test_flow(self):
        try:
            self.mc_client.flush_mc_db()
            self.a2g2_client.flush_a2g2_db()
            job_queue = self.create_job_queue()
            flow_queue = self.create_job_queue()
            self.create_flows()
            self.assertTrue(len(self.mc_client.fetch_tickable_flows()) > 0)
            self.run_flows_to_completion(timeout=60, max_ticks=100,
                                         job_queue=job_queue,
                                         flow_queue=flow_queue)
            self.assert_domain_db_has_expected_state()
        except Exception as exception:
            flows = self.mc_client.fetch_flows()
            flows = [
                {
                    **self._pluck(flow, keys=['created', 'status', 'label']),
                    #'nodes': flow['serialization']['nodes']
                }
                for flow in flows
                if flow['status'] is not 'COMPLETED'
            ]
            self.dump_db_objects(db_objects=flows, label='flow')
            raise

    def create_flows(self):
        flow_spec = reaxys_flow_spec_generator.generate_flow_spec()
        flow_kwargs = {
            'label': flow_spec.get('label'),
            'serialization': FlowEngine.serialize_flow(
                flow=FlowEngine.generate_flow(flow_spec=flow_spec))
        }
        self.mc_client.create_flow(flow_kwargs=flow_kwargs)

    def create_job_queue(self):
        return self.mc_client.create_queue({
            'label': 'reaxys_job_queue',
            'queue_spec': {
                'root_model_spec': 'missions.Job'
            }
        })

    def create_flow_queue(self):
        return self.mc_client.create_queue({
            'label': 'reaxys_job_queue',
            'queue_spec': {
                'root_model_spec': 'missions.Flow'
            }
        })

    def _pluck(self, _dict=None, keys=None):
        keys_set = set(keys)
        return {k:v for k,v in _dict.items() if k in keys_set}

    def assert_domain_db_has_expected_state(self):
        try:
            self.fail("add more tests here")
        except:
            raise

from mc.a2g2.job_engines.a2g2_job_engine import A2G2JobEngine
class FakeJobEngine(object):
    def __init__(self):
        self.wrapped_engine = A2G2JobEngine()

    def execute_command(self, *args, job=None, command=None, **kwargs):
        if job['job_spec']['job_type'].startswith('a2g2.jobs.confgen'):
            pass
        elif job['job_spec']['job_type'].startswith('a2g2.jobs.a2g2_db.query'):
            if command == 'run_job_submission':
                a2g2_db_fixtures = A2G2_DB_Fixtures()
                a2g2_db_fixtures.run_fake_conformer_query(num_results=1)
        else:
            return self.wrapped_engine.execute_command(
                *args, job=job, command=command, **kwargs)

if __name__ == '__main__':
    # Map fake calls to handlers.
    # Will be called in fake execution context, NOT in unittest context.
    parser = argparse.ArgumentParser(description='mocked executables')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    if parsed_args.command == 'fake_qchem':
        QChemFixtures.handle_qchem_command(args=command_args)
    if parsed_args.command == 'fake_confgen':
        ConfgenFixtures.handle_confgen_command(args=command_args)
