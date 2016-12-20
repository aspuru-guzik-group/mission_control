import unittest
from unittest.mock import MagicMock

from ..workflow_engine import WorkflowEngine
from ..workflow import Workflow, BaseNode


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = WorkflowEngine()

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start()
                 for key, patcher in patchers.items()}
        return mocks

class DeserializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_node_classes()

    def setup_node_classes(self):
        class BaseNode(object):
            def __init__(self, *args, id=None, status=None, **kwargs):
                self.id = id
                self.status = status

        self.node_classes = {name: type(name, (BaseNode,), {})
                             for name in ['Type1', 'Type2', 'Type3']}
        for name, node_class in self.node_classes.items():
            self.engine.register_node_class(node_class=node_class)

        self.serialized_workflow = {
            'nodes': [
                {'id': 'a', 'node_type': 'Type1', 'status': 'COMPLETED'},
                {'id': 'b', 'node_type': 'Type2', 'status': 'COMPLETED'},
                {'id': 'c', 'node_type': 'Type1', 'status': 'RUNNING'},
                {'id': 'd', 'node_type' : 'Type2', 'status': 'RUNNING'},
            ],
            'root_node_id': 'a',
            'edges': [
                {'src_id': 'a', 'dest_id': 'b'},
                {'src_id': 'b', 'dest_id': 'c'},
                {'src_id': 'b', 'dest_id': 'd'},
            ]
        }
        self.workflow = self.engine.deserialize_workflow(
            serialized_workflow=self.serialized_workflow)

    def test_has_expected_nodes(self):
        expected_nodes = {}
        for serialized_node in self.serialized_workflow['nodes']:
            expected_node = self.engine.deserialize_node(serialized_node)
            expected_nodes[expected_node.id] = expected_node

        def summarize_nodes(nodes):
            return {
                node.id: {
                    'id': node.id,
                    'node_type': type(node),
                    'status': node.status,
                }
                for node in nodes.values()
            }

        self.assertEqual(summarize_nodes(self.workflow.nodes),
                         summarize_nodes(expected_nodes))

    def test_has_expected_root_node(self):
        self.assertEqual(self.workflow.root_node, self.workflow.nodes['a'])

    def test_has_expected_edges(self):
        expected_edges = {}
        for serialized_edge in self.serialized_workflow['edges']:
            src_id, dest_id = (serialized_edge['src_id'],
                               serialized_edge['dest_id'])
            expected_edges[(src_id, dest_id)] = {
                'src': self.workflow.nodes[src_id],
                'dest': self.workflow.nodes[dest_id]
            }
        self.assertEqual(self.workflow.edges, expected_edges)


class SerializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.workflow = self.generate_workflow()

    def generate_workflow(self):
        workflow = Workflow()
        class BaseNode(object):
            node_type = 'BaseNode'
            def __init__(self, id=None, status=None, state=None, **kwargs):
                self.id = id
                self.status = status
                self.state = state

        nodes = [BaseNode(id=i, status='s_%s' % i, state=i) for i in range(3)]
        workflow.add_nodes(nodes=nodes)
        workflow.root_node = nodes[0]
        edges = [{'src': nodes[0], 'dest': nodes[1]},
                 {'src': nodes[1], 'dest': nodes[2]}]
        workflow.add_edges(edges=edges)
        workflow.state = 'workflow_state'
        workflow.status = 'workflow_status'
        return workflow

    def test_serializes_workflow(self):
        serialization = self.engine.serialize_workflow(self.workflow)
        expected_serialization = {
            'nodes': [{'id': node.id, 'node_type': node.node_type,
                       'status': node.status, 'state': node.state}
                      for node in self.workflow.nodes.values()],
            'root_node_id': self.workflow.root_node.id,
            'edges': [{'src_id': edge['src'].id, 'dest_id': edge['dest'].id}
                      for edge in self.workflow.edges.values()],
            'state': self.workflow.state,
            'status': self.workflow.status,
        }
        self.assertEqual(serialization, expected_serialization)

class TickTestCase(BaseTestCase):
    def test_starts_nearest_pending_nodes(self):
        self.maxDiff = None
        workflow = self.generate_workflow_with_pending_successors()
        precursors = workflow.get_nodes_by_status(status='COMPLETED')
        successors = workflow.get_nodes_by_status(status='PENDING')
        self.assertTrue(len(successors) > 0)

        expected_successor_summaries_before_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 0,
                                                          'status': 'PENDING'})
            for node in successors
        }
        self.assertEqual(self.summarize_nodes(successors),
                         expected_successor_summaries_before_tick)

        self.engine.tick_workflow(workflow)

        expected_successor_summaries_after_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 1,
                                                          'status': 'RUNNING'})
            for node in successors
        }
        self.assertEqual(self.summarize_nodes(successors),
                         expected_successor_summaries_after_tick)
        
        expected_precusor_summaries_after_tick = {
            node.id: self.summarize_node(
                node, overrides={'tick_count': 0, 'status': 'COMPLETED'})
            for node in precursors
        }
        self.assertEqual(self.summarize_nodes(precursors),
                         expected_precusor_summaries_after_tick)

    def summarize_nodes(self, nodes):
        return {node.id: self.summarize_node(node) for node in nodes}

    def summarize_node(self, node, overrides=None):
        if not overrides: overrides = {}
        return {'id': node.id, 'tick_count': node.tick.call_count,
                'status': node.status, **overrides}

    def generate_workflow_with_pending_successors(self):
        workflow = Workflow(root_node=BaseNode(status='COMPLETED'))
        self.add_branch_with_pending_leaf_to_workflow(workflow, depth=3)
        self.add_branch_with_pending_leaf_to_workflow(workflow, depth=2)
        return workflow

    def add_branch_with_pending_leaf_to_workflow(self, workflow, depth=3):
        branch_root = MagicMock()
        branch_root.status = 'COMPLETED'
        workflow.add_child_nodes(parent_node=workflow.root_node,
                                 child_nodes=[branch_root])
        for i in range(depth):
            child_node = MagicMock()
            child_node.status = 'COMPLETED'
            child_node.depth = depth
            workflow.add_child_nodes(parent_node=branch_root,
                                     child_nodes=[child_node])
            branch_root = child_node
        branch_leaf = MagicMock()
        branch_leaf.status = 'PENDING'
        workflow.add_child_nodes(parent_node=branch_root,
                                 child_nodes=[branch_leaf])

    def test_ticks_running_tasks(self):
        workflow = self.generate_workflow_with_running_nodes()
        running_nodes = workflow.get_nodes_by_status(status='RUNNING')
        expected_node_summaries_before_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 0,
                                                          'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_before_tick)
        self.engine.tick_workflow(workflow)
        expected_node_summaries_after_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 1,
                                                          'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_after_tick)

    def generate_workflow_with_running_nodes(self):
        workflow = Workflow(root_node=BaseNode(status='COMPLETED'))
        for i in range(3):
            running_node = MagicMock()
            running_node.status = 'RUNNING'
            workflow.add_node(running_node)
            workflow.connect_nodes(src=workflow.root_node, dest=running_node)
        return workflow

    def test_sets_status_to_completed_if_no_incomplete_nodes(self):
        workflow = Workflow(root_node=BaseNode(status='COMPLETED'))
        self.assertTrue(workflow.status != 'COMPLETED')
        self.engine.tick_workflow(workflow)
        self.assertTrue(workflow.status == 'COMPLETED')

if __name__ == '__main__':
    unittest.main()
