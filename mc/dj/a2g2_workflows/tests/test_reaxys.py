from collections import defaultdict
import unittest
from unittest.mock import patch
from workflow_engines.workflow import Workflow

from ..workflows import reaxys


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

class Reaxys_WorkflowModuleTestCase(BaseTestCase):
    def test_has_expected_node_classes(self):
        self.fail()

class GenerateWorkflowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

        class MockNode(object):
            def __init__(self, *args, **kwargs):
                self.args = args
                self.kwargs = kwargs

            def __eq__(self, other):
                return (self.args == other.args) and \
                        (self.kwargs == other.kwargs)

        self.patchers = {
            'reaxys': patch.multiple(
                reaxys, node_classes=defaultdict(lambda: MockNode))
        }
        self.start_patchers(self.patchers)

    def test_workflow_has_expected_structure(self):
        workflow = reaxys.Reaxys_WorkflowModule.generate_workflow()
        node_classes = reaxys.node_classes
        expected_workflow = Workflow()
        expected_workflow.add_node(
            node_classes['A2G2_DAO'](id='fetch_mol', data=None),
            root=True)
        expected_workflow.add_node(
            node_classes['Reaxys_ConfgenSetupNode'](id='confgen_setup'),
            precursor='fetch_mol')
        expected_workflow.add_node(
            node_classes['Confgen'](id='confgen'),
            precursor='confgen_setup')
        expected_workflow.add_node(
            node_classes['Reaxys_B3LYPSetup'](id='b3lyp_setup'),
            precursor='confgen')
        expected_workflow.add_node(
            node_classes['B3LYPNode'](id='b3lyp'),
            precursor='b3lyp_setup')

        self.assertWorkflowsEqual(workflow, expected_workflow)

    def assertWorkflowsEqual(self, workflow_a, workflow_b):
        self.assertEqual(workflow_a.nodes, workflow_b.nodes)
        self.assertEqual(workflow_a.edges, workflow_b.edges)

if __name__ == '__main__':
    unittest.main()
