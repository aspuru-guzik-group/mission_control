from collections import defaultdict
import unittest
from unittest.mock import call, DEFAULT, patch, MagicMock

from ..workflow import Workflow


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
            nodes.append(self.generate_node(id=i))
        return nodes

    def generate_node(self, **kwargs):
        node = MagicMock()
        node.configure_mock(**kwargs)
        return node

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
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()

    def test_adds_node(self):
        self.assertFalse(self.node.id in self.workflow.nodes)
        self.workflow.add_node(node=self.node)
        self.assertEqual(self.workflow.nodes[self.node.id], self.node)

    def test_sets_workflow_on_node(self):
        self.workflow.add_node(node=self.node)
        self.assertEqual(self.node.workflow, self.workflow)

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
    def setUp(self):
        super().setUp()
        self.nodes = self.generate_nodes(n=2)
        self.edge = {'src': self.nodes[0], 'dest': self.nodes[1]}
        self.expected_edge_key = (self.edge['src'].id, self.edge['dest'].id)

    def test_add_edge(self):
        self.assertFalse(self.expected_edge_key in self.workflow.edges)
        self.workflow.add_edge(edge=self.edge)
        self.assertEqual(self.workflow.edges[self.expected_edge_key], self.edge)

    def test_updates_edges_by_node_id(self):
        self.assertTrue(self.nodes[0].id not in self.workflow.edges_by_node_id)
        self.assertTrue(self.nodes[1].id not in self.workflow.edges_by_node_id)
        self.workflow.add_edge(edge=self.edge)
        for node in self.nodes:
            self.assertEqual(
                self.workflow.edges_by_node_id[node.id][self.expected_edge_key],
                self.edge)

class ConnectNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'workflow': patch.multiple(self.workflow, add_edge=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_connect_nodes(self):
        nodes = self.generate_nodes(n=2)
        self.workflow.add_nodes(nodes=nodes)
        self.workflow.connect_nodes(src=nodes[0], dest=nodes[1])
        self.assertEqual(self.mocks['workflow']['add_edge'].call_args,
                         call(edge={'src': nodes[0], 'dest': nodes[1]}))

class GetNearestPendingNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_workflow()

    def setup_workflow(self):
        self.workflow.root_node = self.generate_node(id='ROOT',
                                                     status='COMPLETED')
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
            self.add_linear_branch_to_workflow(workflow=self.workflow,
                                               branch_nodes=branch_nodes)

    def add_linear_branch_to_workflow(self, workflow=None, branch_nodes=None):
        tail = workflow.root_node
        for node in branch_nodes:
            workflow.add_nodes(nodes=[node])
            workflow.connect_nodes(src=tail, dest=node)
            tail = node

    def test_get_nearest_pending_nodes(self):
        expected_nearest_pending_nodes = [self.linear_branches['1'][1],
                                          self.linear_branches['2'][2]]
        self.assertEqual(set(self.workflow.get_nearest_pending_nodes()),
                         set(expected_nearest_pending_nodes))

class GetChildNodesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_workflow(workflow=self.workflow)

    def setup_workflow(self, workflow=None):
        workflow.root_node = self.generate_node(id='ROOT')
        node_1 = self.generate_node(id='1')
        self.workflow.add_node(node_1)
        workflow.connect_nodes(src=workflow.root_node, dest=node_1)
        node_1_1 = self.generate_node(id='1_1')
        self.workflow.add_node(node_1_1)
        workflow.connect_nodes(src=node_1, dest=node_1_1)
        node_1_2 = self.generate_node(id='1_2')
        self.workflow.add_node(node_1_2)
        workflow.connect_nodes(src=node_1, dest=node_1_2)

    def test_gets_child_nodes(self):
        expected_child_nodes = [self.workflow.nodes['1_1'],
                                self.workflow.nodes['1_2']]
        child_nodes = self.workflow.get_child_nodes(
            parent_node=self.workflow.nodes['1'])
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
            self.workflow.add_node(node)
        _filter = lambda node: node.answer == 'yes'
        result = self.workflow.filter_nodes(filters=[_filter])
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
                self.workflow.add_node(node)
        for status in statuses:
            result = self.workflow.get_nodes_by_status(status=status)
            expected_result = nodes[status]
            self.assertEqual(set(result), set(expected_result))

class HasIncompleteNodesTestCase(BaseTestCase):
    def test_has_incomplete(self):
        self.workflow.add_node(self.generate_node(status='PENDING'))
        self.assertTrue(self.workflow.has_incomplete_nodes())

    def test_does_not_have_incomplete(self):
        self.workflow.add_node(self.generate_node(status='COMPLETED'))
        self.assertFalse(self.workflow.has_incomplete_nodes())


if __name__ == '__main__':
    unittest.main()
