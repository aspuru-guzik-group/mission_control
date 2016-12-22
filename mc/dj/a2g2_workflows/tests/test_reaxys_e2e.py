from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from workflow_engines.workflow_engine import WorkflowEngine

from ..workflow.reaxys import generate_reaxys_workflow


class ReaxysWorkflowTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx = self.generate_ctx()
        self.engine = self.generate_workflow_engine(ctx=self.ctx)
        self.mol_key = 'some mol key'
        self.workflow = self.generate_workflow(mol_key=self.mol_key)

    def generate_ctx(self):
        return defaultdict(MagicMock)

    def generate_workflow_engine(self, ctx=None):
        engine = WorkflowEngine(ctx=ctx)
        return engine

    def generate_workflow(self, mol=None):
        return generate_reaxys_workflow(mol=mol)

    def test_run_workflow(self):
        states = []
        states.append(self.get_state(prev_state=None))

        def tick_and_capture_state(label=None):
            print("tick_and_capture", label)
            self.engine.tick_workflow(workflow=self.workflow)
            states.append(self.get_state(prev_state=states[-1]))

        # Set initial data.
        root_node = self.workflow.root_node
        self.assertEqual(root_node.data['output']['mol_key'], self.mol_key)
        self.assertEqual(root_node.status, 'COMPLETED')

        # Fetch the mol.
        fetch_mol_node = self.workflow.nodes['fetch_mol']
        self.self.assertEqual(fetch_mol_node.status, 'PENDING')
        tick_and_capture_state()
        self.assertEqual(self.ctx['a2g2_dao'].Mol.objects.get.call_args,
                         call(key=self.mol_key))
        self.assertEaual(fetch_mol_node.status, 'COMPLETED')
        expected_mol = self.ctx['a2g2_dao'].Mol.objects.get.return_value
        self.assertEaual(fetch_mol_node.data['output']['mol'], expected_mol)

        # Wire confgen.
        wire_confgen_node = self.workflow.nodes['wire_confgen']
        self.assertEaual(wire_confgen_node.status, 'PENDING')
        tick_and_capture_state()
        self.assertEaual(wire_confgen_node.status, 'COMPLETED')
        confgen_node = self.workflow.nodes['confgen']
        self.assertEaual(confgen_node.status, 'PENDING')
        self.assertEaual(confgen_node.data['input']['smiles'],
                         fetch_mol_node.data['output']['mol']['smiles'])

        # Run confgen.
        tick_and_capture_state()
        self.assertEaual(confgen_node.status, 'RUNNING')
        self.assertEqual(self.ctx['create_job'].call_args, {
            'type': 'confgen',
            'spec': {
                'smiles': confgen_node.data['input']['smiles']
            }
        })
        expected_confgen_job_id = self.ctx['create_job'].return_value
        self.assertEqual(confgen_node.data['_job_id'], expected_confgen_job_id)

        # Complete confgen job.
        completed_confgen_job = {'status': 'COMPLETED',
                                 'data': {'output': {'output_dir': 'some dir'}}}
        self.workflow.jobs[expected_confgen_job_id] = completed_confgen_job
        tick_and_capture_state()
        self.assertEqual(confgen_node.status, 'COMPLETED')
        self.assertEqual(confgen_node.data['output']['output'],
                         completed_confgen_job['data']['output'])

        # Wire ingest_confgen.
        tick_and_capture_state()
        wire_ingest_confgen_node = self.workflow.nodes['wire_ingest_confgen']
        self.assertEqual(wire_ingest_confgen_node.status, 'COMPLETED')
        # ingest_confgen
        ingest_confgen_node = self.workflow.nodes['ingest_confgen']
        self.assertEqual(ingest_confgen_node.status, 'PENDING')
        self.assertEqual(ingest_confgen_node.data['input']['input_dir'],
                         confgen_node['data']['output']['output_dir'])
        tick_and_capture_state()

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
