import unittest
from unittest.mock import call, MagicMock
from uuid import uuid4

from flow_engines.flow_engine import FlowEngine

from ..flows.reaxys import Reaxys_Flow_Module


class ReaxysFlowTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx = self.generate_ctx()
        self.setup_create_job_side_effect()
        self.engine = self.generate_flow_engine(ctx=self.ctx)
        self.mol_key = 'some mol key'
        self.flow = self.generate_flow(mol_key=self.mol_key)

    def generate_ctx(self):
        return {
            'a2g2_dao': MagicMock(),
            'create_job': MagicMock(),
        }

    def setup_create_job_side_effect(self):
        """For tracking job ids, for faking job completion."""
        self.job_ids = []
        def create_job_side_effect(*args, **kwargs):
            job_id = uuid4()
            self.job_ids.append(job_id)
            return job_id
        self.ctx['create_job'].side_effect = create_job_side_effect

    def generate_flow_engine(self, ctx=None):
        engine = FlowEngine(ctx=ctx)
        return engine

    def generate_flow(self, mol_key=None):
        return Reaxys_Flow_Module.generate_flow(mol_key=mol_key)

    def test_run_flow(self):

        # Helpers
        self.tick_counter = 0
        def _tick():
            print("tick_counter: %s" % self.tick_counter)
            self.tick_counter += 1
            self.engine.tick_flow(flow=self.flow, ctx=self.ctx)
        nodes = self.flow.nodes

        self.ctx['a2g2_dao'].query.return_value = {
            'mols': [{'smiles': 'some smiles'}]
        }

        # t = 0
        self.assertEqual(nodes['fetch_mol'].data['input'],
                         {'query': {
                             'object': 'Mol',
                             'filters': [('key', self.mol_key)]
                         }})
        self.assertEqual(nodes['fetch_mol'].status, 'PENDING')

        # t = 1
        _tick()
        self.assertEqual(nodes['fetch_mol'].status, 'COMPLETED')
        self.assertEqual(nodes['fetch_mol'].data['output']['result_set'],
                         self.ctx['a2g2_dao'].query.return_value)
        self.assertEqual(nodes['confgen_setup'].status, 'PENDING')

        # t = 2
        _tick()
        self.assertEqual(nodes['confgen_setup'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen'].status, 'PENDING')
        self.assertEqual(nodes['confgen'].data['input']['smiles'],
                         nodes['fetch_mol'].data['output']['result_set']\
                         ['mols'][0]['smiles'])

        # t = 3
        _tick()
        self.assertEqual(nodes['confgen'].status, 'RUNNING')
        self.assertEqual(
            self.ctx['create_job'].call_args,
            call(job_kwargs={
                'type': 'confgen',
                'spec': {
                    'smiles': nodes['confgen'].data['input']['smiles'],
                    'confgen': nodes['confgen'].data['input']['confgen_spec'],
                }
            })
        )
        confgen_job_id = self.job_ids[-1]
        self.assertEqual(nodes['confgen'].job_node.data['job_id'],
                         confgen_job_id)

        # Fake confgen job completion.
        confgen_job = {'status': 'COMPLETED',
                                 'data': {'output': {'dir': 'some dir'}}}
        self.flow.jobs[confgen_job_id] = confgen_job

        # t = 4
        _tick()
        self.assertEqual(nodes['confgen'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen'].data['output'],
                         confgen_job['data']['output'])
        self.assertEqual(nodes['confgen_parse_setup'].status, 'PENDING')

        # t = 5
        _tick()
        self.assertEqual(nodes['confgen_parse_setup'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen_parse'].status, 'PENDING')
        self.assertEqual(
            nodes['confgen_parse'].data['input'],
            {
                'parse': {
                    'dir': nodes['confgen'].data['output']['dir'],
                }
            }
        )

        # t = 6
        _tick()
        confgen_parse_job_id = self.job_ids[-1]
        self.assertEqual(nodes['confgen_parse'].job_node.data['job_id'],
                         confgen_parse_job_id)

        # Fake confgen_parse job completion.
        confgen_parse_job = {'status': 'COMPLETED',
                             'data': {'output': {'dir': 'some dir'}}}
        self.flow.jobs[confgen_parse_job_id] = confgen_parse_job

        # t = 7
        _tick()
        self.assertEqual(nodes['confgen_parse'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen_parse'].data['output'],
                         confgen_parse_job['data']['output'])
        self.assertEqual(nodes['confgen_ingest_setup'].status, 'PENDING')

        # t = 8
        _tick()
        self.assertEqual(nodes['confgen_ingest_setup'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen_ingest'].status, 'PENDING')
        self.assertEqual(
            nodes['confgen_ingest'].data['input'],
            {
                'ingest': {
                    'object': 'Geom',
                    'item_kwargs': {'mol': self.mol_key},
                    'dir': nodes['confgen_parse'].data['output']['dir'],
                }
            }
        )
        self.fail()

        # t = 9
        _tick()

        # B3LYP Task
        # Fetches geom from db 
        # Runs b3lyp node.
        # Generates b3lyp job.
        # <fake completion of b3lyp job>
        # Ingests confgen results.
        # <fake the ingest>

        # Flow is completed.

    def get_state(self, prev_state=None):
        state = {}
        return state

if __name__ == '__main__':
    unittest.main()
