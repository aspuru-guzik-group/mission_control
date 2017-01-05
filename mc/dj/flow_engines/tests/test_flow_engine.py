import unittest
from unittest.mock import Mock

from ..flow_engine import FlowEngine
from ..flow import Flow, BaseNode


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = FlowEngine()

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
            def __init__(self, *args, id=None, status=None, ctx=None, **kwargs):
                self.id = id
                self.status = status
                self.ctx = ctx

        self.node_classes = {name: type(name, (BaseNode,), {})
                             for name in ['Type1', 'Type2', 'Type3']}
        for name, node_class in self.node_classes.items():
            self.engine.register_node_class(node_class=node_class)

        self.serialized_flow = {
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
            ],
            'status': 'status',
        }
        self.flow = self.engine.deserialize_flow(
            serialized_flow=self.serialized_flow)

    def test_has_expected_nodes(self):
        expected_nodes = {}
        for serialized_node in self.serialized_flow['nodes']:
            expected_node = self.engine.deserialize_node(serialized_node)
            expected_nodes[expected_node.id] = expected_node

        def summarize_nodes(nodes):
            return {
                node.id: {
                    'id': node.id,
                    'node_type': type(node),
                    'status': node.status,
                    'ctx': node.ctx,
                }
                for node in nodes.values()
            }

        self.assertEqual(summarize_nodes(self.flow.nodes),
                         summarize_nodes(expected_nodes))

    def test_has_expected_root_node(self):
        self.assertEqual(self.flow.root_node, self.flow.nodes['a'])

    def test_has_expected_edges(self):
        expected_edges = {}
        for serialized_edge in self.serialized_flow['edges']:
            src_id, dest_id = (serialized_edge['src_id'],
                               serialized_edge['dest_id'])
            expected_edges[(src_id, dest_id)] = {
                'src': self.flow.nodes[src_id],
                'dest': self.flow.nodes[dest_id]
            }
        self.assertEqual(self.flow.edges, expected_edges)

    def test_has_expected_status(self):
        self.assertEqual(
            self.flow.status, self.serialized_flow['status'])


class SerializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = self.generate_flow()

    def generate_flow(self):
        flow = Flow()
        class BaseNode(object):
            node_type = 'BaseNode'
            def __init__(self, id=None, status=None, state=None, **kwargs):
                self.id = id
                self.status = status
                self.state = state

        nodes = [BaseNode(id=i, status='s_%s' % i, state=i) for i in range(3)]
        flow.add_nodes(nodes=nodes)
        flow.root_node = nodes[0]
        edges = [{'src': nodes[0], 'dest': nodes[1]},
                 {'src': nodes[1], 'dest': nodes[2]}]
        flow.add_edges(edges=edges)
        flow.state = 'flow_state'
        flow.status = 'flow_status'
        return flow

    def test_serializes_flow(self):
        serialization = self.engine.serialize_flow(self.flow)
        expected_serialization = {
            'nodes': [{'id': node.id, 'node_type': node.node_type,
                       'status': node.status, 'state': node.state}
                      for node in self.flow.nodes.values()],
            'root_node_id': self.flow.root_node.id,
            'edges': [{'src_id': edge['src'].id, 'dest_id': edge['dest'].id}
                      for edge in self.flow.edges.values()],
            'state': self.flow.state,
            'status': self.flow.status,
        }
        self.assertEqual(serialization, expected_serialization)

class TickTestCase(BaseTestCase):
    def test_starts_nearest_pending_nodes(self):
        self.maxDiff = None
        flow = self.generate_flow_with_pending_successors()
        precursors = flow.get_nodes_by_status(status='COMPLETED')
        successors = flow.get_nodes_by_status(status='PENDING')
        self.assertTrue(len(successors) > 0)

        expected_successor_summaries_before_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 0,
                                                          'status': 'PENDING'})
            for node in successors
        }
        self.assertEqual(self.summarize_nodes(successors),
                         expected_successor_summaries_before_tick)

        self.engine.tick_flow(flow)

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

    def generate_flow_with_pending_successors(self):
        flow = Flow()
        root_node = Mock()
        root_node.configure_mock(id='ROOT', status='COMPLETED')
        flow.add_node(root_node, as_root=True)
        self.add_branch_with_pending_leaf_to_flow(flow, depth=3)
        self.add_branch_with_pending_leaf_to_flow(flow, depth=2)
        return flow

    def add_branch_with_pending_leaf_to_flow(self, flow, depth=3):
        branch_root = Mock()
        branch_root.status = 'COMPLETED'
        flow.add_node(node=branch_root, precursor=flow.root_node)
        for i in range(depth):
            child_node = Mock()
            child_node.status = 'COMPLETED'
            child_node.depth = depth
            flow.add_node(node=child_node, precursor=branch_root)
            branch_root = child_node
        branch_leaf = Mock()
        branch_leaf.status = 'PENDING'
        flow.add_node(node=branch_leaf, precursor=branch_root)

    def test_ticks_running_tasks(self):
        flow = self.generate_flow_with_running_nodes()
        running_nodes = flow.get_nodes_by_status(status='RUNNING')
        expected_node_summaries_before_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 0,
                                                          'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_before_tick)
        self.engine.tick_flow(flow)
        expected_node_summaries_after_tick = {
            node.id: self.summarize_node(node, overrides={'tick_count': 1,
                                                          'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_after_tick)

    def generate_flow_with_running_nodes(self):
        flow = Flow()
        flow.add_node(BaseNode(status='COMPLETED'), as_root=True)
        for i in range(3):
            running_node = Mock()
            running_node.configure_mock(id=i, status='RUNNING')
            flow.add_node(running_node, precursor=flow.root_node)
        return flow

    def test_sets_status_to_completed_if_no_incomplete_nodes(self):
        flow = Flow(root_node=BaseNode(status='COMPLETED'))
        self.assertTrue(flow.status != 'COMPLETED')
        self.engine.tick_flow(flow)
        self.assertTrue(flow.status == 'COMPLETED')

if __name__ == '__main__':
    unittest.main()
