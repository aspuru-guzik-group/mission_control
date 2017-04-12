import unittest
from unittest.mock import call, DEFAULT, patch, MagicMock
from uuid import uuid4

from .. import flow_engine
from ..flow import Flow


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task_handler = MagicMock
        self.engine = flow_engine.FlowEngine(task_handler=self.task_handler)

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
        self.setup_flow()

    def setup_flow(self):
        self.serialized_flow = {
            'data': 'some data',
            'nodes': [
                {'node_key': 'a', 'status': 'COMPLETED'},
                {'node_key': 'b', 'status': 'COMPLETED'},
                {'node_key': 'c', 'status': 'RUNNING'},
                {'node_key': 'd', 'status': 'RUNNING'},
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

    def test_has_expected_nodes(self):
        self.assertEqual(
            self.flow.nodes,
            {node['node_key']: node for node in self.serialized_flow['nodes']}
        )

    def test_has_expected_root_node_key(self):
        self.assertEqual(self.flow.root_node_key,
                         self.flow.nodes['a']['node_key'])

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
        flow.data = 'some data'
        nodes = [{'key': i, 'status': 's_%s' % i} for i in range(3)]
        for node in nodes: flow.add_node(node=node)
        flow.root_node_key = nodes[0]['node_key']
        edges = [{'src_key': nodes[0]['node_key'], 'dest_key': nodes[1]['key']},
                 {'src_key': nodes[1]['node_key'], 'dest_key': nodes[2]['key']}]
        for edge in edges: flow.add_edge(edge=edge)
        flow.status = 'flow_status'
        return flow

    def test_serializes_basic_attrs(self):
        for attr in ['data', 'root_node_key', 'status']:
            self.assertEqual(self.serialization[attr], getattr(self.flow, attr))

    def test_serializes_nodes(self):
        def sorted_node_list(node_list):
            return sorted(node_list, key=lambda node: node['node_key'])
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
        def mock_tick_node(*args, flow=None, node=None, flow_ctx=None):
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
        flow.add_node(node={'node_key': 'ROOT', 'status': 'COMPLETED'},
                      as_root=True)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=3)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=2)
        return flow

    def _add_branch_with_pending_leaf_to_flow(self, flow, depth=3):
        branch_root = flow.add_node(node={'status': 'COMPLETED'},
                                    precursor_keys=[flow.root_node_key])
        for i in range(depth):
            child_node = flow.add_node(
                node={'node_key': str(uuid4()), 'status': 'COMPLETED',
                      'depth': depth},
                precursor_keys=[branch_root['node_key']]
            )
            branch_root = child_node
        flow.add_node(node={'status': 'PENDING'},
                      precursor_keys=[branch_root['node_key']])

    def test_ticks_running_nodes(self):
        self.maxDiff = None
        flow = self._generate_flow_with_running_nodes()
        running_nodes = flow.get_nodes_by_status(status='RUNNING')
        expected_node_summaries_before_tick = {
            node['node_key']: self.summarize_node(
                node, overrides={'tick_count': None, 'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_before_tick)
        self.engine.tick_flow(flow)
        expected_node_summaries_after_tick = {
            node['node_key']: self.summarize_node(node, overrides={'tick_count': 1,
                                                           'status': 'RUNNING'})
            for node in running_nodes
        }
        self.assertEqual(self.summarize_nodes(running_nodes),
                         expected_node_summaries_after_tick)

    def _generate_flow_with_running_nodes(self):
        flow = Flow()
        flow.add_node(node={'status': 'COMPLETED'}, as_root=True)
        for i in range(3):
            flow.add_node(node={'node_key': i, 'status': 'RUNNING'},
                          precursor_keys=[flow.root_node_key])
        return flow

    def summarize_nodes(self, nodes):
        return {node['node_key']: self.summarize_node(node) for node in nodes}

    def summarize_node(self, node, overrides=None):
        if not overrides: overrides = {}
        return {
            **{attr: node.get(attr, None)
               for attr in ['node_key', 'tick_count', 'status']},
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
        self.flow_ctx = MagicMock()
        self.engine.get_task_runner_for_node = MagicMock()
        self.expected_task_runner = \
                self.engine.get_task_runner_for_node.return_value

    def test_calls_task_runner(self):
        self.engine.tick_node(flow=self.flow, node=self.node,
                              flow_ctx=self.flow_ctx)
        self.assertEqual(
            len(self.expected_task_runner.tick_tasks.call_args_list), 1)

    def test_calls_complete_node_if_node_completes(self):
        self.engine.complete_node = MagicMock()
        self.expected_task_runner.tick_tasks.return_value = 'COMPLETED'
        self.engine.tick_node(flow=self.flow, node=self.node, flow_ctx=self.flow_ctx)
        self.assertEqual(self.engine.complete_node.call_args,
                         call(flow=self.flow, node=self.node))

class GetTaskRunnerForNodeTestCase(BaseTestCase):
    @patch.object(flow_engine, 'BaseTaskRunner')
    def test_generates_task_runner_with_expected_kwargs(self, MockTaskRunner):
        node = MagicMock()
        flow = MagicMock()
        flow_ctx = MagicMock()
        self.engine.get_task_runner_for_node(node=node, flow=flow,
                                             flow_ctx=flow_ctx)
        call_kwargs = MockTaskRunner.call_args[1]
        self.assertEqual(call_kwargs['get_tasks'](), node.get('node_tasks'))
        self.assertEqual(
            call_kwargs['get_task_context'](),
            self.engine.get_task_context(node=node, flow=flow,
                                         flow_ctx=flow_ctx)
        )
        self.assertEqual(call_kwargs['task_handler'], self.engine.task_handler)

class CompleteNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.node = {}

    def test_sets_node_status_to_completed(self):
        self.engine.complete_node(flow=self.flow, node=self.node)
        self.assertEqual(self.node['status'], 'COMPLETED')

class CompleteFlowTestCase(BaseTestCase):
    def test_sets_status_to_completed(self):
        flow = Flow()
        self.assertTrue(flow.status != 'COMPLETED')
        self.engine.complete_flow(flow=flow)
        self.assertTrue(flow.status == 'COMPLETED')

class StartNodeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.node = {}

    def test_sets_status_to_running(self):
        self.engine.start_node(flow=self.flow, node=self.node)
        self.assertEqual(self.node['status'], 'RUNNING')

if __name__ == '__main__':
    unittest.main()
