import unittest
from unittest.mock import call, DEFAULT, patch, MagicMock, Mock
from uuid import uuid4

from ..flow_engine import FlowEngine
from ..flow import Flow


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
        self.setup_node_engines()

    def setup_node_engines(self):
        class MockNodeEngine(object):
            def __init__(self, *args, **kwargs):
                self.__dict__.update(kwargs)

        self.node_engines = {name: type(name, (MockNodeEngine,), {})
                             for name in ['Engine1', 'Engine2', 'Engine3']}
        for name, node_engine in self.node_engines.items():
            self.engine.register_node_engine(node_engine=node_engine)

        self.serialized_flow = {
            'data': 'some data',
            'input': 'some input',
            'output': 'some output',
            'nodes': [
                {'key': 'a', 'node_engine': 'Engine1', 'status': 'COMPLETED'},
                {'key': 'b', 'node_engine': 'Engine2', 'status': 'COMPLETED'},
                {'key': 'c', 'node_engine': 'Engine1', 'status': 'RUNNING'},
                {'key': 'd', 'node_engine': 'Engine2', 'status': 'RUNNING'},
            ],
            'root_node_key': 'a',
            'edges': [
                {'src_key': 'a', 'dest_key': 'b'},
                {'src_key': 'b', 'dest_key': 'c'},
                {'src_key': 'b', 'dest_key': 'd'},
            ],
            'status': 'status',
        }
        self.flow = self.engine.deserialize_flow(
            serialized_flow=self.serialized_flow)

    def test_has_expected_data(self):
        self.assertEqual(self.flow.data, self.serialized_flow['data'])

    def test_has_expected_input(self):
        self.assertEqual(self.flow.input, self.serialized_flow['input'])

    def test_has_expected_output(self):
        self.assertEqual(self.flow.output, self.serialized_flow['output'])

    def test_has_expected_nodes(self):
        self.assertEqual(
            self.flow.nodes,
            {node['key']: node for node in self.serialized_flow['nodes']}
        )

    def test_has_expected_root_node_key(self):
        self.assertEqual(self.flow.root_node_key, self.flow.nodes['a']['key'])

    def test_has_expected_edges(self):
        expected_edges = {}
        for edge in self.serialized_flow['edges']:
            edge_key = (edge['src_key'], edge['dest_key'])
            expected_edges[edge_key] = edge
        self.assertEqual(self.flow.edges, expected_edges)

    def test_has_expected_status(self):
        self.assertEqual(
            self.flow.status, self.serialized_flow['status'])

class SerializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = self.generate_flow()
        self.serialization = self.engine.serialize_flow(self.flow)

    def generate_flow(self):
        flow = Flow()

        flow.input = 'some input'
        flow.output = 'some output'
        flow.data = 'some data'

        nodes = [{'key': i, 'status': 's_%s' % i} for i in range(3)]
        for node in nodes: flow.add_node(node=node)
        flow.root_node_key = nodes[0]['key']
        edges = [{'src_key': nodes[0]['key'], 'dest_key': nodes[1]['key']},
                 {'src_key': nodes[1]['key'], 'dest_key': nodes[2]['key']}]
        for edge in edges: flow.add_edge(edge=edge)
        flow.status = 'flow_status'
        return flow

    def test_serializes_basic_attrs(self):
        for attr in ['data', 'input', 'output', 'root_node_key', 'status']:
            self.assertEqual(self.serialization[attr], getattr(self.flow, attr))

    def test_serializes_nodes(self):
        def sorted_node_list(node_list):
            return sorted(node_list, key=lambda node: node['key'])
        self.assertEqual(sorted_node_list(self.serialization['nodes']),
                         sorted_node_list(self.flow.nodes.values()))

    def test_serializes_edges(self):
        def sorted_edge_list(edge_list):
            def get_edge_key_str(edge):
                return "%s,%s" % (edge['src_key'], edge['dest_key'])
            return sorted(edge_list, key=get_edge_key_str)
        self.assertEqual(sorted_edge_list(self.serialization['edges']),
                         sorted_edge_list(self.flow.edges.values()))

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {'engine': patch.multiple(
            self.engine, start_node=DEFAULT, complete_flow=DEFAULT)}
        self.mocks = self.start_patchers(patchers=self.patchers)
        def mock_tick_node(*args, flow=None, node=None, ctx=None):
            node.setdefault('tick_count', 0)
            node['tick_count'] += 1
        self.engine.tick_node = mock_tick_node

    def test_starts_pending_flow(self):
        flow = Flow(status='PENDING')
        self.engine.tick_flow(flow=flow)
        self.assertEqual(flow.status, 'RUNNING')

    def test_starts_nearest_pending_nodes(self):
        flow = self._generate_flow_with_pending_successors()
        successors = flow.get_nodes_by_status(status='PENDING')
        self.assertTrue(len(successors) > 0)
        self.engine.tick_flow(flow)
        expected_call_args_list = [call(flow=flow, node=successor)
                                   for successor in successors]
        self.assertEqual(
            sorted(self.mocks['engine']['start_node'].call_args_list,
                   key=lambda c: id(c[1]['node'])),
            sorted(expected_call_args_list,
                   key=lambda c: id(c[2]['node'])))

    def _generate_flow_with_pending_successors(self):
        flow = Flow()
        flow.add_node(key='ROOT', node={'status': 'COMPLETED'},
                      as_root=True)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=3)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=2)
        return flow

    def _add_branch_with_pending_leaf_to_flow(self, flow, depth=3):
        branch_root = flow.add_node(node={'status': 'COMPLETED'},
                                    precursor_keys=[flow.root_node_key])
        for i in range(depth):
            child_node = flow.add_node(
                key=str(uuid4()),
                node={'status': 'COMPLETED', 'depth': depth},
                precursor_keys=[branch_root['key']]
            )
            branch_root = child_node
        flow.add_node(node={'status': 'PENDING'},
                      precursor_keys=[branch_root['key']])

    def test_ticks_running_nodes(self):
        self.maxDiff = None
        flow = self._generate_flow_with_running_nodes()
        running_nodes = flow.get_nodes_by_status(status='RUNNING')
        expected_node_summaries_before_tick = {
            node['key']: self.summarize_node(
                node, overrides={'tick_count': None, 'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_before_tick)
        self.engine.tick_flow(flow)
        expected_node_summaries_after_tick = {
            node['key']: self.summarize_node(node, overrides={'tick_count': 1,
                                                           'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_after_tick)

    def _generate_flow_with_running_nodes(self):
        flow = Flow()
        flow.add_node(node={'status': 'COMPLETED'}, as_root=True)
        for i in range(3):
            flow.add_node(key=i, node={'status': 'RUNNING'},
                          precursor_keys=[flow.root_node_key])
        return flow

    def summarize_nodes(self, nodes):
        return {node['key']: self.summarize_node(node) for node in nodes}

    def summarize_node(self, node, overrides=None):
        if not overrides: overrides = {}
        return {
            **{attr: node.get(attr, None)
               for attr in ['key', 'tick_count', 'status']},
            **overrides
        }

    def test_calls_complete_flow_if_no_incomplete_nodes(self):
        flow = Flow()
        flow.add_node(node={'status': 'COMPLETED'}, as_root=True)
        self.assertEqual(self.mocks['engine']['complete_flow'].call_count, 0)
        self.engine.tick_flow(flow)
        self.assertEqual(self.mocks['engine']['complete_flow'].call_count, 1)

class TickNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.node = self.flow.add_node(node={})
        self.engine.process_actions = Mock()
        self.engine.get_engine_for_node = Mock()
        self.mock_engine = self.engine.get_engine_for_node.return_value
        self.ctx = MagicMock()
        self.engine.complete_node = Mock()

    def test_calls_node_engine_tick_node(self):
        self.engine.tick_node(flow=self.flow, node=self.node, ctx=self.ctx)
        self.assertEqual(
            self.mock_engine.tick_node.call_args,
            call(node=self.node, ctx={**self.ctx, 'flow': self.flow}))

    def test_calls_complete_node_if_node_completes(self):
        def mock_node_completion(*args, **kwargs):
            self.node['status'] = 'COMPLETED'
        self.mock_engine.tick_node.side_effect = mock_node_completion
        self.engine.tick_node(flow=self.flow, node=self.node, ctx=self.ctx)
        self.assertEqual(self.engine.complete_node.call_args,
                         call(flow=self.flow, node=self.node))

class CompleteNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.node = self.flow.add_node(node={
            'post_complete_actions': ['action_%s' % i for i in range(3)]
        })
        self.engine.process_actions = Mock()

    def test_processes_post_complete_actions(self):
        self.engine.complete_node(flow=self.flow, node=self.node)
        self.assertEqual(self.engine.process_actions.call_args,
                         call(actions=self.node['post_complete_actions'],
                              ctx={'node': self.node, 'flow': self.flow}))

class CompleteFlowTestCase(BaseTestCase):
    def test_sets_status_to_completed(self):
        flow = Flow()
        self.assertTrue(flow.status != 'COMPLETED')
        self.engine.complete_flow(flow=flow)
        self.assertTrue(flow.status == 'COMPLETED')

    def test_sets_output_for_single_tail_node(self):
        flow = Flow()
        tail_node = flow.add_node(node={'status': 'COMPLETED'},
                                  as_root=True)
        tail_node['output'] = 'some output'
        self.engine.complete_flow(flow=flow)
        self.assertEqual(flow.output, tail_node['output'])

    def test_sets_output_for_multiple_tail_nodes(self):
        flow = Flow()
        root_node = flow.add_node(node={'status': 'COMPLETED'},
                                  as_root=True)
        tail_nodes = []
        for i in range(3):
            tail_node = flow.add_node(node={'status': 'COMPLETED'},
                                      precursor_keys=[root_node['key']])
            tail_node['output'] = 'output_%s' % i
            tail_nodes.append(tail_node)
        self.engine.complete_flow(flow=flow)
        expected_output = {tail_node['key']: tail_node['output']
                           for tail_node in tail_nodes}
        self.assertEqual(flow.output, expected_output)

class StartNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.node = self.flow.add_node(node={
            'pre_start_actions': ['action_%s' % i for i in range(3)]
        })
        self.engine.process_actions = Mock()

    def test_processes_pre_start_actions(self):
        self.engine.start_node(flow=self.flow, node=self.node)
        self.assertEqual(self.engine.process_actions.call_args,
                         call(actions=self.node['pre_start_actions'],
                              ctx={'node': self.node, 'flow': self.flow}))

    def test_sets_status_to_running(self):
        self.engine.start_node(flow=self.flow, node=self.node)
        self.assertEqual(self.node['status'], 'RUNNING')

class ProcessActionsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.actions = ['action_%s' % i for i in range(3)]
        self.engine.action_processor = Mock()
        self.ctx = Mock()

    def test_calls_action_processor_for_each_action(self):
        self.engine.process_actions(actions=self.actions, ctx=self.ctx)
        self.assertEqual(
            self.engine.action_processor.process_action.call_args_list,
            [call(action=action, ctx=self.ctx) for action in self.actions]
        )

class GenerateFlowTestCase(BaseTestCase):
    def test_generates_flow_from_spec(self):
        class MockFlowGenerator(object):
            flow_type = 'some flow type'
            generate_flow = Mock()
        self.engine.register_flow_generator_class(MockFlowGenerator)
        flow_spec = {'flow_type': MockFlowGenerator.flow_type}
        flow = self.engine.generate_flow(flow_spec=flow_spec)
        self.assertEqual(MockFlowGenerator.generate_flow.call_args,
                         call(flow_spec=flow_spec))
        self.assertEqual(flow, MockFlowGenerator.generate_flow.return_value)

if __name__ == '__main__':
    unittest.main()
