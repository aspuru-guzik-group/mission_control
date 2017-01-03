import unittest
from workflow_engines.workflow import Workflow

from ..workflows import reaxys
from ..nodes.a2g2_dao import A2G2_DAO_Node
from ..nodes.b3lyp import B3LYP_Node
from ..nodes.confgen import Confgen_Node


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()

    def tearDown(self):
        super().tearDown()
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def start_patchers(self, patchers):
        return {key: patcher.start() for key, patcher in patchers.items()}

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

class Reaxys_Workflow_ModuleTestCase(BaseTestCase):
    def test_has_expected_node_classes(self):
        node_classes = reaxys.Reaxys_Workflow_Module.get_node_classes()
        expected_node_classes = {
            'A2G2_DAO': A2G2_DAO_Node,
            'Reaxys_Confgen_Setup': reaxys.Reaxys_Confgen_Setup_Node,
            'Confgen': Confgen_Node,
            'Reaxys_B3LYP_Setup': reaxys.Reaxys_B3LYP_Setup_Node,
            'B3LYP': B3LYP_Node,
        }
        self.assertEqual(node_classes, expected_node_classes)

class GenerateWorkflowTestCase(BaseTestCase):
    def test_workflow_has_expected_structure(self):
        workflow = reaxys.Reaxys_Workflow_Module.generate_workflow()
        node_classes = reaxys.Reaxys_Workflow_Module.get_node_classes()
        expected_workflow = Workflow()
        expected_workflow.add_node(
            node=node_classes['A2G2_DAO'](id='fetch_mol', data=None),
            as_root=True)
        expected_workflow.add_node(
            node=node_classes['Reaxys_Confgen_Setup'](id='confgen_setup'),
            precursor='fetch_mol')
        expected_workflow.add_node(
            node=node_classes['Confgen'](id='confgen'),
            precursor='confgen_setup')
        expected_workflow.add_node(
            node=node_classes['Reaxys_B3LYP_Setup'](id='b3lyp_setup'),
            precursor='confgen')
        expected_workflow.add_node(
            node=node_classes['B3LYP'](id='b3lyp'),
            precursor='b3lyp_setup')

        self.assertWorkflowsEqual(workflow, expected_workflow)

    def assertWorkflowsEqual(self, workflow_a, workflow_b):
        self.assertEqual(workflow_a.nodes.keys(), workflow_b.nodes.keys())
        self.assertEqual(workflow_a.edges.keys(), workflow_b.edges.keys())

if __name__ == '__main__':
    unittest.main()
