import unittest
from flow_engines.flow import Flow

from ..flows.rpi import RPI_Flow_Module


class RunParseIngest_TestCase(unittest.TestCase):
    def test_flow_has_expected_structure(self):
        flow = rpi.RPI_Flow_Module.generate_flow()
        node_classes = rpi.RPI_Flow_Module.get_node_classes()
        expected_flow = Flow()
        expected_flow.add_node(
            node=node_classes['run'](id='run'),
            as_root=True)
        expected_flow.add_node(
            node=node_classes['rpi:parse_setup'](id='parse_setup'),
            precursor='run')
        expected_flow.add_node(
            node=node_classes['rpi:parse'](id='parse'),
            precursor='parse_setup')
        expected_flow.add_node(
            node=node_classes['rpi:ingest_setup'](id='ingest_setup'),
            precursor='parse')
        expected_flow.add_node(
            node=node_classes['rpi:ingest'](id='ingest'),
            precursor='ingest_setup')
        self.assertFlowsEqual(flow, expected_flow)

    def assertFlowsEqual(self, flow_a, flow_b):
        self.assertEqual(flow_a.nodes.keys(), flow_b.nodes.keys())
        self.assertEqual(flow_a.edges.keys(), flow_b.edges.keys())

if __name__ == '__main__':
    unittest.main()
