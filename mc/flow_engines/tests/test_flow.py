from collections import defaultdict
import unittest
from unittest.mock import Mock
from uuid import uuid4

from ..flow import Flow


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow = Flow()

    def generate_nodes(self, n=3):
        nodes = []
        for i in range(n):
            nodes.append(self.generate_node(node_key=i))
        return nodes

    def generate_node(self, node_key=None, **kwargs):
        if node_key is None: node_key = str(uuid4())
        node = {'node_key': node_key, **kwargs}
        return node

    def _node_list_to_dict(self, node_list):
        return {node['node_key']: node for node in node_list}

class AddNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node(node_key='node')

    def test_adds_node(self):
        self.assertFalse(self.node['node_key'] in self.flow.nodes)
        self.flow.add_node(node=self.node)
        self.assertEqual(self.flow.nodes[self.node['node_key']], self.node)

    def test_returns_node(self):
        result = self.flow.add_node(node=self.node)
        self.assertEqual(result, self.node)

    def test_generates_key_if_unset(self):
        self.flow.generate_node_key = Mock()
        del self.node['node_key']
        self.flow.add_node(node=self.node)
        self.assertEqual(self.node['node_key'],
                         self.flow.generate_node_key.return_value)

    def test_adds_edges_for_precursors(self):
        precursors = [self.flow.add_node(node=self.generate_node(key=i))
                      for i in range(1)]
        self.flow.add_node(
            node=self.node,
            precursor_keys=[precursor['node_key'] for precursor in precursors])
        for precursor in precursors:
            self.assertTrue(self.flow.has_edge(src_key=precursor['node_key'],
                                               dest_key=self.node['node_key']))

    def test_adds_edges_for_successor(self):
        successors = [self.flow.add_node(node=self.generate_node(key=i))
                      for i in range(3)]
        self.flow.add_node(
            node=self.node,
            successor_keys=[successor['node_key'] for successor in successors])
        for successor in successors:
            self.assertTrue(self.flow.has_edge(src_key=self.node['node_key'],
                                               dest_key=successor['node_key']))

class AddEdgeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.nodes = [self.flow.add_node(node=node)
                      for node in self.generate_nodes(n=2)]
        self.edge = {'src_key': self.nodes[0]['node_key'],
                     'dest_key': self.nodes[1]['node_key']}
        self.expected_edge_key = (self.edge['src_key'], self.edge['dest_key'])

    def test_add_edge(self):
        self.assertFalse(self.expected_edge_key in self.flow.edges)
        self.flow.add_edge(edge=self.edge)
        self.assertEqual(self.flow.edges[self.expected_edge_key], self.edge)

    def test_updates_edges_by_node_key(self):
        self.assertTrue(self.nodes[0]['node_key'] not in self.flow.edges_by_node_key)
        self.assertTrue(self.nodes[1]['node_key'] not in self.flow.edges_by_node_key)
        self.flow.add_edge(edge=self.edge)
        self.assertEqual(
            self.flow.edges_by_node_key[self.edge['src_key']]['outgoing'],
            {self.expected_edge_key: self.edge})
        self.assertEqual(
            self.flow.edges_by_node_key[self.edge['dest_key']]['incoming'],
            {self.expected_edge_key: self.edge})

class GetNearestPendingNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow()

    def setup_flow(self):
        self.linear_branches = {
            '1': [
                self.generate_node(node_key='1.1', status='COMPLETED'),
                self.generate_node(node_key='1.2', status='PENDING'),
                self.generate_node(node_key='1.3', status='PENDING')
            ],
            '2': [
                self.generate_node(node_key='2.1', status='COMPLETED'),
                self.generate_node(node_key='2.2', status='COMPLETED'),
                self.generate_node(node_key='2.3', status='PENDING'),
                self.generate_node(node_key='2.4', status='PENDING'),
            ],
        }
        for branch_nodes in self.linear_branches.values():
            self.add_linear_branch_to_flow(flow=self.flow,
                                           branch_nodes=branch_nodes)

    def add_linear_branch_to_flow(self, flow=None, branch_nodes=None):
        tail = flow.nodes[flow.ROOT_NODE_KEY]
        for node in branch_nodes:
            tail = flow.add_node(node=node, precursor_keys=[tail['node_key']])

    def test_get_nearest_pending_nodes(self):
        expected_nearest_pending_nodes = [self.linear_branches['1'][1],
                                          self.linear_branches['2'][2]]
        self.assertEqual(self.flow.get_nearest_pending_nodes(),
                         expected_nearest_pending_nodes)

class GetSuccessorsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow(flow=self.flow)

    def setup_flow(self, flow=None):
        node_1 = self.generate_node(node_key='1')
        self.flow.add_node(node=node_1, precursor_keys=[flow.ROOT_NODE_KEY])
        node_1_1 = self.generate_node(node_key='1_1')
        self.flow.add_node(node=node_1_1, precursor_keys=[node_1['node_key']])
        node_1_2 = self.generate_node(node_key='1_2')
        self.flow.add_node(node=node_1_2, precursor_keys=[node_1['node_key']])

    def test_gets_successors(self):
        successors = self.flow.get_successors(node=self.flow.nodes['1'])
        expected_successors = [self.flow.nodes['1_1'], self.flow.nodes['1_2']]
        self.assertEqual(self._node_list_to_dict(successors),
                         self._node_list_to_dict(expected_successors))

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
            self.flow.add_node(node=node)
        _filter = lambda node: node['answer'] == 'yes'
        result = self.flow.filter_nodes(filters=[_filter])
        expected_result = nodes['yes']
        self.assertEqual(self._node_list_to_dict(result),
                         self._node_list_to_dict(expected_result))

class GetNodesByStatusTestCase(BaseTestCase):
    def test_gets_nodes_by_status(self):
        nodes = defaultdict(list)
        statuses = ['status_%s' % i for i in range(3)]
        for status in statuses:
            for i in range(3):
                node = self.generate_node(status=status)
                nodes[status].append(node)
                self.flow.add_node(node=node)
        for status in statuses:
            result = self.flow.get_nodes_by_status(status=status)
            expected_result = nodes[status]
            self.assertEqual(self._node_list_to_dict(result),
                             self._node_list_to_dict(expected_result))

class HasIncompleteNodesTestCase(BaseTestCase):
    def test_has_incomplete(self):
        self.flow.add_node(node=self.generate_node(status='PENDING'))
        self.assertTrue(self.flow.has_incomplete_nodes())

    def test_does_not_have_incomplete(self):
        self.flow.add_node(node=self.generate_node(status='COMPLETED'))
        self.assertFalse(self.flow.has_incomplete_nodes())

class GetTailNodesTestCase(BaseTestCase):
    def test_gets_tail_nodes(self):
        tail_nodes = [
            self.flow.add_node(node=self.generate_node(),
                               precursor_keys=[self.flow.ROOT_NODE_KEY])
            for i in range(3)]
        self.assertTrue(self.flow.get_tail_nodes(), tail_nodes)

    def test_handles_flow_with_only_root_node(self):
        self.assertTrue(self.flow.get_tail_nodes(),
                        self.flow.nodes[self.flow.ROOT_NODE_KEY])

if __name__ == '__main__':
    unittest.main()
