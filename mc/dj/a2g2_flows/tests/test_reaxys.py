import unittest
from flow_engines.flow import Flow

from ..flows import reaxys
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

class Reaxys_Flow_ModuleTestCase(BaseTestCase):
    def test_has_expected_node_classes(self):
        node_classes = reaxys.Reaxys_Flow_Module.get_node_classes()
        expected_node_classes = {
            'A2G2_DAO': A2G2_DAO_Node,
            'Reaxys_Confgen_Setup': reaxys.Reaxys_Confgen_Setup_Node,
            'Confgen': Confgen_Node,
            'Reaxys_B3LYP_Setup': reaxys.Reaxys_B3LYP_Setup_Node,
            'B3LYP': B3LYP_Node,
        }
        self.assertEqual(node_classes, expected_node_classes)

class GenerateFlowTestCase(BaseTestCase):
    def test_flow_has_expected_structure(self):
        flow = reaxys.Reaxys_Flow_Module.generate_flow()
        node_classes = reaxys.Reaxys_Flow_Module.get_node_classes()
        expected_flow = Flow()
        expected_flow.add_node(
            node=node_classes['A2G2_DAO'](id='fetch_mol', data=None),
            as_root=True)
        expected_flow.add_node(
            node=node_classes['Reaxys_Confgen_Setup'](id='confgen_setup'),
            precursor='fetch_mol')
        expected_flow.add_node(
            node=node_classes['Confgen'](id='confgen'),
            precursor='confgen_setup')
        expected_flow.add_node(
            node=node_classes['Reaxys_B3LYP_Setup'](id='b3lyp_setup'),
            precursor='confgen')
        expected_flow.add_node(
            node=node_classes['B3LYP'](id='b3lyp'),
            precursor='b3lyp_setup')

        self.assertFlowsEqual(flow, expected_flow)

    def assertFlowsEqual(self, flow_a, flow_b):
        self.assertEqual(flow_a.nodes.keys(), flow_b.nodes.keys())
        self.assertEqual(flow_a.edges.keys(), flow_b.edges.keys())

if __name__ == '__main__':
    unittest.main()
