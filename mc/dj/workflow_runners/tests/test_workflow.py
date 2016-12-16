import unittest
from unittest.mock import call, DEFAULT, patch

from ..workflow import Workflow

class StubNode(object):
    def __init__(self, id=None):
        self.id = id


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.workflow = Workflow()

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start() for key, patcher in patchers.items()}
        return mocks

    def generate_nodes(self, n=3):
        nodes = []
        for i in range(n):
            nodes.append(StubNode(id=i))
        return nodes

class AddNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'workflow': patch.multiple(self.workflow, add_node=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_adds_nodes(self):
        nodes = [i for i in range(3)]
        self.workflow.add_nodes(nodes=nodes)
        self.assertEqual(self.mocks['workflow']['add_node'].call_args_list,
                         [call(node=node) for node in nodes])

class AddNodeTestCase(BaseTestCase):
    def test_adds_node(self):
        node = StubNode(id='id')
        self.assertFalse(node.id in self.workflow.nodes)
        self.workflow.add_node(node=node)
        self.assertEqual(self.workflow.nodes[node.id], node)

class AddEdgesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'workflow': patch.multiple(self.workflow, add_edge=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_adds_edges(self):
        edges = [i for i in range(3)]
        self.workflow.add_edges(edges=edges)
        self.assertEqual(self.mocks['workflow']['add_edge'].call_args_list,
                         [call(edge=edge) for edge in edges])

class AddEdgeTestCase(BaseTestCase):
    def test_add_edge(self):
        nodes = self.generate_nodes(n=2)
        edge = {'src': nodes[0], 'dest': nodes[1]}
        expected_edge_key = (edge['src'].id, edge['dest'].id)
        self.assertFalse(expected_edge_key in self.workflow.edges)
        self.workflow.add_edge(edge=edge)
        self.assertEqual(self.workflow.edges[expected_edge_key], edge)

if __name__ == '__main__':
    unittest.main()
