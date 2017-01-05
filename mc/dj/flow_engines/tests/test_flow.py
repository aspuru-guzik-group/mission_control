from collections import defaultdict
import unittest
from unittest.mock import call, DEFAULT, patch, Mock
from uuid import uuid4

from ..flow import Flow


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow = Flow()

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
            nodes.append(self.generate_node(id=i))
        return nodes

    def generate_node(self, id=None, **kwargs):
        if id is None: id = str(uuid4())
        node = Mock()
        node.configure_mock(id=id, **kwargs)
        return node

class AddNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'flow': patch.multiple(self.flow, add_node=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_adds_nodes(self):
        nodes = self.generate_nodes(n=3)
        self.flow.add_nodes(nodes=nodes)
        self.assertEqual(self.mocks['flow']['add_node'].call_args_list,
                         [call(node=node) for node in nodes])

class AddNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node(id='node')

    def test_adds_node(self):
        self.assertFalse(self.node.id in self.flow.nodes)
        self.flow.add_node(node=self.node)
        self.assertEqual(self.flow.nodes[self.node.id], self.node)

    def test_returns_node(self):
        result = self.flow.add_node(node=self.node)
        self.assertEqual(result, self.node)

    def test_sets_flow_on_node(self):
        self.flow.add_node(node=self.node)
        self.assertEqual(self.node.flow, self.flow)

    def test_add_as_root(self):
        self.flow.add_node(node=self.node, as_root=True)
        self.assertEqual(self.node.flow, self.flow)
        self.assertEqual(self.flow.root_node, self.node)

    def test_add_with_single_precursor(self):
        precursor = self.flow.add_node(node=self.generate_node())
        self.flow.add_node(node=self.node, precursor=precursor)
        self.assertTrue(self.flow.has_edge(src=precursor, dest=self.node))

    def test_add_with_precursor_id(self):
        precursor = self.generate_node(id='precursor_id')
        self.flow.add_node(node=precursor)
        self.flow.add_node(node=self.node, precursor='precursor_id')
        self.assertTrue(self.flow.has_edge(src=precursor, dest=self.node))

    def test_add_with_many_precursors(self):
        precursors = [self.flow.add_node(node=self.generate_node(id=i))
                      for i in range(1)]
        self.flow.add_node(node=self.node, precursor=precursors)
        for precursor in precursors:
            self.assertTrue(
                self.flow.has_edge(src=precursor, dest=self.node))

    def test_add_with_single_successor(self):
        successor = self.flow.add_node(node=self.generate_node())
        self.flow.add_node(node=self.node, successor=successor)
        self.assertTrue(self.flow.has_edge(src=self.node, dest=successor))

    def test_add_with_successsor_id(self):
        successor = self.generate_node(id='successor_id')
        self.flow.add_node(node=successor)
        self.flow.add_node(node=self.node, successor='successor_id')
        self.assertTrue(self.flow.has_edge(src=self.node, dest=successor))

    def test_add_with_many_successors(self):
        successors = [self.flow.add_node(node=self.generate_node(id=i))
                      for i in range(3)]
        self.flow.add_node(node=self.node, successor=successors)
        for successor in successors:
            self.assertTrue(
                self.flow.has_edge(src=self.node, dest=successor))

class AddEdgesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'flow': patch.multiple(self.flow, add_edge=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_adds_edges(self):
        edges = [{'src': i, 'dest': i} for i in range(3)]
        self.flow.add_edges(edges=edges)
        self.assertEqual(self.mocks['flow']['add_edge'].call_args_list,
                         [call(**edge) for edge in edges])

class AddEdgeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.nodes = self.generate_nodes(n=2)
        self.flow.add_nodes(nodes=self.nodes)
        self.edge = {'src': self.nodes[0], 'dest': self.nodes[1]}
        self.expected_edge_key = (self.edge['src'].id, self.edge['dest'].id)

    def test_add_edge(self):
        self.assertFalse(self.expected_edge_key in self.flow.edges)
        self.flow.add_edge(**self.edge)
        self.assertEqual(self.flow.edges[self.expected_edge_key], self.edge)

    def test_updates_edges_by_node_id(self):
        self.assertTrue(self.nodes[0].id not in self.flow.edges_by_node_id)
        self.assertTrue(self.nodes[1].id not in self.flow.edges_by_node_id)
        self.flow.add_edge(**self.edge)
        for node in self.nodes:
            self.assertEqual(
                self.flow.edges_by_node_id[node.id][self.expected_edge_key],
                self.edge)

class GetNearestPendingNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow()

    def setup_flow(self):
        self.flow.add_node(
            self.generate_node(id='ROOT', status='COMPLETED'), as_root=True)
        self.linear_branches = {
            '1': [
                self.generate_node(id='1.1', status='COMPLETED'),
                self.generate_node(id='1.2', status='PENDING'),
                self.generate_node(id='1.3', status='PENDING')
            ],
            '2': [
                self.generate_node(id='2.1', status='COMPLETED'),
                self.generate_node(id='2.2', status='COMPLETED'),
                self.generate_node(id='2.3', status='PENDING'),
                self.generate_node(id='2.4', status='PENDING'),
            ],
        }
        for branch_nodes in self.linear_branches.values():
            self.add_linear_branch_to_flow(flow=self.flow,
                                               branch_nodes=branch_nodes)

    def add_linear_branch_to_flow(self, flow=None, branch_nodes=None):
        tail = flow.root_node
        for node in branch_nodes:
            tail = flow.add_node(node=node, precursor=[tail])

    def test_get_nearest_pending_nodes(self):
        expected_nearest_pending_nodes = [self.linear_branches['1'][1],
                                          self.linear_branches['2'][2]]
        self.assertEqual(set(self.flow.get_nearest_pending_nodes()),
                         set(expected_nearest_pending_nodes))

class GetChildNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow(flow=self.flow)

    def setup_flow(self, flow=None):
        self.flow.add_node(self.generate_node(id='ROOT'), as_root=True)
        node_1 = self.generate_node(id='1')
        self.flow.add_node(node_1, precursor=flow.root_node)
        node_1_1 = self.generate_node(id='1_1')
        self.flow.add_node(node_1_1, precursor=node_1)
        node_1_2 = self.generate_node(id='1_2')
        self.flow.add_node(node_1_2, precursor=node_1)

    def test_gets_child_nodes(self):
        expected_child_nodes = [self.flow.nodes['1_1'],
                                self.flow.nodes['1_2']]
        child_nodes = self.flow.get_child_nodes(
            parent_node=self.flow.nodes['1'])
        self.assertEqual(set(child_nodes), set(expected_child_nodes))

class FilterNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_filters(self):
        nodes = defaultdict(list)
        for i in range(4):
            if i % 2: answer = 'yes'
            else: answer = 'no'
            node = self.generate_node(answer=answer)
            nodes[answer].append(node)
            self.flow.add_node(node)
        _filter = lambda node: node.answer == 'yes'
        result = self.flow.filter_nodes(filters=[_filter])
        expected_result = nodes['yes']
        self.assertEqual(set(result), set(expected_result))

class GetNodesByStatusTestCase(BaseTestCase):
    def test_gets_nodes_by_status(self):
        nodes = defaultdict(list)
        statuses = ['status_%s' % i for i in range(3)]
        for status in statuses:
            for i in range(3):
                node = self.generate_node(status=status)
                nodes[status].append(node)
                self.flow.add_node(node)
        for status in statuses:
            result = self.flow.get_nodes_by_status(status=status)
            expected_result = nodes[status]
            self.assertEqual(set(result), set(expected_result))

class HasIncompleteNodesTestCase(BaseTestCase):
    def test_has_incomplete(self):
        self.flow.add_node(self.generate_node(status='PENDING'))
        self.assertTrue(self.flow.has_incomplete_nodes())

    def test_does_not_have_incomplete(self):
        self.flow.add_node(self.generate_node(status='COMPLETED'))
        self.assertFalse(self.flow.has_incomplete_nodes())


if __name__ == '__main__':
    unittest.main()
