import unittest
from unittest.mock import call, MagicMock

from workflow_engines.workflow_engine import WorkflowEngine

from ..workflows.reaxys import Reaxys_Workflow_Module


class ReaxysWorkflowTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx = self.generate_ctx()
        self.engine = self.generate_workflow_engine(ctx=self.ctx)
        self.mol_key = 'some mol key'
        self.workflow = self.generate_workflow(mol_key=self.mol_key)

    def generate_ctx(self):
        return {
            'a2g2_dao': MagicMock(),
            'create_job': MagicMock(),
        }

    def generate_workflow_engine(self, ctx=None):
        engine = WorkflowEngine(ctx=ctx)
        return engine

    def generate_workflow(self, mol_key=None):
        return Reaxys_Workflow_Module.generate_workflow(mol_key=mol_key)

    def test_run_workflow(self):

        # Helpers
        self.tick_counter = 0
        def _tick():
            print("tick_counter: %s" % self.tick_counter)
            self.tick_counter += 1
            self.engine.tick_workflow(workflow=self.workflow, ctx=self.ctx)
        nodes = self.workflow.nodes

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
        expected_confgen_job_id = self.ctx['create_job'].return_value
        self.assertEqual(nodes['confgen'].job_node.data['job_id'],
                         expected_confgen_job_id)

        # Fake confgen job completion.
        completed_confgen_job = {'status': 'COMPLETED',
                                 'data': {'output': {'output_dir': 'some dir'}}}
        self.workflow.jobs[expected_confgen_job_id] = completed_confgen_job

        # t = 4
        _tick()
        self.assertEqual(nodes['confgen'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen'].data['output'],
                         completed_confgen_job['data']['output'])

        # t = 5
        _tick()
        self.assertEqual(nodes['confgen_ingest_setup'].status, 'COMPLETED')
        self.assertEqual(nodes['confgen_ingest'].status, 'PENDING')
        self.assertEqual(
            nodes['confgen_ingest'].data['input'],
            {
                'ingest': {
                    'object': 'Geom',
                    'item_kwargs': {'mol': self.mol_key},
                    'xyz_dir': nodes['confgen'].data['output']['output_dir'],
                }
            }
        )
        self.fail()

        # t = 6
        _tick()

        # B3LYP Task
        # Fetches geom from db 
        # Runs b3lyp node.
        # Generates b3lyp job.
        # <fake completion of b3lyp job>
        # Ingests confgen results.
        # <fake the ingest>

        # Workflow is completed.

    def get_state(self, prev_state=None):
        state = {}
        return state

if __name__ == '__main__':
    unittest.main()
