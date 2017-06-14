from collections import defaultdict
import json
import unittest
from unittest.mock import call, DEFAULT, patch, MagicMock
from uuid import uuid4

from .. import flow_engine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task_handler = MagicMock()
        self.engine = flow_engine.FlowEngine(task_handler=self.task_handler)

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start()
                 for key, patcher in patchers.items()}
        return mocks

class GenerateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_spec = defaultdict(MagicMock, **{
            'data': MagicMock(),
            'task_spec': [MagicMock() for i in range(3)],
        })

    @patch.object(flow_engine, 'Flow')
    def test_generates_flow_from_flow_spec(self, MockFlow):
        flow = self.engine.generate_flow(flow_spec=self.flow_spec)
        self.assertEqual(flow, MockFlow.return_value)
        self.assertEqual(MockFlow.call_args,
                         call(data=self.flow_spec['data']))
        expected_add_task_call_args_list = [
            call(**task_spec) for task_spec in self.flow_spec['task_specs']
        ]
        self.assertEqual(flow.add_task.call_args_list,
                         expected_add_task_call_args_list)

class DeserializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_dict = {
            'data': 'some data',
            'tasks': {
                'a': {'key': 'a', 'status': 'COMPLETED'},
                'b': {'key': 'b', 'status': 'COMPLETED'},
                'c': {'key': 'c', 'status': 'RUNNING'},
                'd': {'key': 'd', 'status': 'RUNNING'},
            },
            'edges': [
                {'src_key': 'ROOT', 'dest_key': 'a'},
                {'src_key': 'a', 'dest_key': 'b'},
                {'src_key': 'b', 'dest_key': 'c'},
                {'src_key': 'b', 'dest_key': 'd'},
            ],
            'status': 'status',
        }

    def _deserialize_flow(self):
        return self.engine.deserialize_flow(
            serialized_flow=json.dumps(self.flow_dict))

    def test_has_expected_data(self):
        flow = self._deserialize_flow()
        self.assertEqual(flow.data, self.flow_dict['data'])

    def test_has_expected_tasks(self):
        flow = self._deserialize_flow()
        expected_tasks = {task['key']: task
                          for task in self.flow_dict['tasks'].values()}
        expected_tasks[flow.ROOT_TASK_KEY] = flow.tasks[flow.ROOT_TASK_KEY]
        self.assertEqual(set(flow.tasks.keys()), set(expected_tasks.keys()))

    def test_has_expected_edges(self):
        flow = self._deserialize_flow()
        expected_edges = {}
        for edge in self.flow_dict['edges']:
            edge_key = (edge['src_key'], edge['dest_key'])
            expected_edges[edge_key] = edge
        self.assertEqual(flow.edges, expected_edges)

    def test_has_expected_status(self):
        flow = self._deserialize_flow()
        self.assertEqual(flow.status, self.flow_dict['status'])

    def test_data_defaults_to_dict(self):
        self.flow_dict = {}
        flow = self._deserialize_flow()
        self.assertEqual(flow.data, {})

class SerializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = self.generate_flow()
        self.flow_dict = json.loads(self.engine.serialize_flow(self.flow))

    def generate_flow(self):
        flow = flow_engine.Flow()
        flow.data = 'some data'
        tasks = [{'key': 'key_%s' % i, 'status': 's_%s' % i} for i in range(3)]
        for task in tasks: flow.add_task(task=task)
        edges = [{'src_key': tasks[0]['key'], 'dest_key': tasks[1]['key']},
                 {'src_key': tasks[1]['key'], 'dest_key': tasks[2]['key']}]
        for edge in edges: flow.add_edge(edge=edge)
        flow.status = 'flow_status'
        return flow

    def test_serializes_basic_attrs(self):
        for attr in flow_engine.FlowEngine.simple_flow_serialization_attrs:
            self.assertEqual(self.flow_dict[attr], getattr(self.flow, attr))

    def test_serializes_tasks(self):
        expected_serialized_tasks = {
            key: task for key, task in self.flow.tasks.items()
            if key != self.flow.ROOT_TASK_KEY
        }
        self.assertEqual(self.flow_dict['tasks'], expected_serialized_tasks)

    def test_serializes_edges(self):
        def sorted_edge_list(edge_list):
            def get_edge_key_str(edge):
                return "%s,%s" % (edge['src_key'], edge['dest_key'])
            return sorted(edge_list, key=get_edge_key_str)
        self.assertEqual(sorted_edge_list(self.flow_dict['edges']),
                         sorted_edge_list(self.flow.edges.values()))

class RunFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = MagicMock()
        self.engine.tick_flow = MagicMock()
        self.tick_counter = 0

    def _setup_mock_tick(self, flow_end_status=None, num_ticks=3, max_ticks=5):
        def mock_tick(*args, flow=None, **kwargs):
            self.tick_counter += 1
            if self.tick_counter == num_ticks: flow.status = flow_end_status
            if self.tick_counter > max_ticks:
                raise Exception("exceeded max ticks")
        self.engine.tick_flow.side_effect = mock_tick

    def test_stops_if_flow_completes(self):
        self._setup_mock_tick(flow_end_status='COMPLETED')
        self.engine.run_flow(flow=self.flow)
        self.assertEqual(self.flow.status, 'COMPLETED')

    def test_stops_if_flow_fails(self):
        self._setup_mock_tick(flow_end_status='FAILED')
        self.engine.run_flow(flow=self.flow)
        self.assertEqual(self.flow.status, 'FAILED')

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {'engine': patch.multiple(
            self.engine, start_task=DEFAULT, complete_flow=DEFAULT)}
        self.mocks = self.start_patchers(patchers=self.patchers)
        def mock_tick_task(*args, flow=None, task=None, task_context=None):
            task.setdefault('tick_count', 0)
            task['tick_count'] += 1
        self.engine.tick_task = mock_tick_task

    def test_starts_pending_flow(self):
        flow = flow_engine.Flow(status='PENDING')
        self.engine.tick_flow(flow=flow)
        self.assertEqual(flow.status, 'RUNNING')

    def test_starts_nearest_pending_tasks(self):
        flow = self._generate_flow_with_pending_successors()
        successors = flow.get_tasks_by_status(status='PENDING')
        self.assertTrue(len(successors) > 0)
        self.engine.tick_flow(flow)
        expected_call_args_list = [call(flow=flow, task=successor)
                                   for successor in successors]
        self.assertEqual(
            sorted(self.mocks['engine']['start_task'].call_args_list,
                   key=lambda c: id(c[1]['task'])),
            sorted(expected_call_args_list,
                   key=lambda c: id(c[2]['task'])))

    def _generate_flow_with_pending_successors(self):
        flow = flow_engine.Flow()
        self._add_branch_with_pending_leaf_to_flow(flow, depth=3)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=2)
        return flow

    def _add_branch_with_pending_leaf_to_flow(self, flow, depth=3):
        branch_root = flow.add_task(task={'status': 'COMPLETED',
                                          'precursors': [flow.ROOT_TASK_KEY]})
        for i in range(depth):
            child_task = flow.add_task(
                task={'key': str(uuid4()), 'status': 'COMPLETED',
                      'precursors': [branch_root['key']], 'depth': depth}
            )
            branch_root = child_task
            flow.add_task(task={'status': 'PENDING',
                                'precursors': [branch_root['key']]})

    def test_ticks_running_tasks(self):
        self.maxDiff = None
        flow = self._generate_flow_with_running_tasks()
        running_tasks = flow.get_tasks_by_status(status='RUNNING')
        expected_task_summaries_before_tick = {
            task['key']: self.summarize_task(
                task, overrides={'tick_count': None, 'status': 'RUNNING'})
            for task in running_tasks
        }
        self.assertEqual(self.summarize_tasks(running_tasks),
                         expected_task_summaries_before_tick)
        self.engine.tick_flow(flow)
        expected_task_summaries_after_tick = {
            task['key']: self.summarize_task(task, overrides={'tick_count': 1,
                                                           'status': 'RUNNING'})
            for task in running_tasks
        }
        self.assertEqual(self.summarize_tasks(running_tasks),
                         expected_task_summaries_after_tick)

    def _generate_flow_with_running_tasks(self):
        flow = flow_engine.Flow()
        for i in range(3):
            flow.add_task(task={'key': i, 'status': 'RUNNING',
                                'precursors': [flow.ROOT_TASK_KEY]})
        return flow

    def summarize_tasks(self, tasks):
        return {task['key']: self.summarize_task(task) for task in tasks}

    def summarize_task(self, task, overrides=None):
        if not overrides: overrides = {}
        return {
            **{attr: task.get(attr, None)
               for attr in ['key', 'tick_count', 'status']},
            **overrides
        }

    def test_calls_complete_flow_if_no_incomplete_tasks(self):
        flow = flow_engine.Flow()
        self.assertEqual(self.mocks['engine']['complete_flow'].call_count, 0)
        self.engine.tick_flow(flow)
        self.assertEqual(self.mocks['engine']['complete_flow'].call_count, 1)

class TickTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = flow_engine.Flow()
        self.task = self.flow.add_task(task={})
        self.task_context = MagicMock()

    def test_calls_complete_task_if_task_completes(self):
        self.engine.complete_task = MagicMock()
        def set_task_status_to_completed(task, *args, **kwargs):
            task['status'] = 'COMPLETED'
        self.engine.task_handler.tick_task.side_effect = \
                set_task_status_to_completed
        self.engine.tick_task(flow=self.flow, task=self.task,
                              task_context=self.task_context)
        self.assertEqual(self.engine.complete_task.call_args,
                         call(flow=self.flow, task=self.task))

class CompleteTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = flow_engine.Flow()
        self.task = {}

    def test_sets_task_status_to_completed(self):
        self.engine.complete_task(flow=self.flow, task=self.task)
        self.assertEqual(self.task['status'], 'COMPLETED')

class CompleteFlowTestCase(BaseTestCase):
    def test_sets_status_to_completed_if_no_errors(self):
        flow = flow_engine.Flow()
        self.engine.complete_flow(flow=flow)
        self.assertEqual(flow.status, 'COMPLETED')

    def test_fails_if_has_errors(self):
        flow = flow_engine.Flow()
        flow.data['errors'] = MagicMock()
        self.engine.complete_flow(flow=flow)
        self.assertEqual(flow.status, 'FAILED')

class StartTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = flow_engine.Flow()
        self.task = {}

    def test_sets_status_to_running(self):
        self.engine.start_task(flow=self.flow, task=self.task)
        self.assertEqual(self.task['status'], 'RUNNING')

if __name__ == '__main__':
    unittest.main()
