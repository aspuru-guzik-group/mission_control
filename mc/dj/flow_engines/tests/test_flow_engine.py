import unittest
from unittest.mock import call, DEFAULT, patch, Mock
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
        self.setup_task_engines()

    def setup_task_engines(self):
        class MockTaskEngine(object):
            def __init__(self, *args, **kwargs):
                self.__dict__.update(kwargs)

        self.task_engines = {name: type(name, (MockTaskEngine,), {})
                             for name in ['Engine1', 'Engine2', 'Engine3']}
        for name, task_engine in self.task_engines.items():
            self.engine.register_task_engine(task_engine=task_engine)

        self.serialized_flow = {
            'data': 'some data',
            'input': 'some input',
            'output': 'some output',
            'tasks': [
                {'key': 'a', 'task_engine': 'Engine1', 'status': 'COMPLETED'},
                {'key': 'b', 'task_engine': 'Engine2', 'status': 'COMPLETED'},
                {'key': 'c', 'task_engine': 'Engine1', 'status': 'RUNNING'},
                {'key': 'd', 'task_engine': 'Engine2', 'status': 'RUNNING'},
            ],
            'root_task_key': 'a',
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

    def test_has_expected_tasks(self):
        self.assertEqual(
            self.flow.tasks,
            {task['key']: task for task in self.serialized_flow['tasks']}
        )

    def test_has_expected_root_task_key(self):
        self.assertEqual(self.flow.root_task_key, self.flow.tasks['a']['key'])

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

        tasks = [{'key': i, 'status': 's_%s' % i} for i in range(3)]
        for task in tasks: flow.add_task(task=task)
        flow.root_task_key = tasks[0]['key']
        edges = [{'src_key': tasks[0]['key'], 'dest_key': tasks[1]['key']},
                 {'src_key': tasks[1]['key'], 'dest_key': tasks[2]['key']}]
        for edge in edges: flow.add_edge(edge=edge)
        flow.status = 'flow_status'
        return flow

    def test_serializes_basic_attrs(self):
        for attr in ['data', 'input', 'output', 'root_task_key', 'status']:
            self.assertEqual(self.serialization[attr], getattr(self.flow, attr))

    def test_serializes_tasks(self):
        def sorted_task_list(task_list):
            return sorted(task_list, key=lambda task: task['key'])
        self.assertEqual(sorted_task_list(self.serialization['tasks']),
                         sorted_task_list(self.flow.tasks.values()))

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
            self.engine, start_task=DEFAULT, complete_flow=DEFAULT)}
        self.mocks = self.start_patchers(patchers=self.patchers)
        def mock_tick_task(*args, task=None, ctx=None):
            task.setdefault('tick_count', 0)
            task['tick_count'] += 1
        self.engine.tick_task = mock_tick_task

    def test_starts_pending_flow(self):
        flow = Flow(status='PENDING')
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
        flow = Flow()
        flow.add_task(key='ROOT', task={'status': 'COMPLETED'},
                      as_root=True)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=3)
        self._add_branch_with_pending_leaf_to_flow(flow, depth=2)
        return flow

    def _add_branch_with_pending_leaf_to_flow(self, flow, depth=3):
        branch_root = flow.add_task(task={'status': 'COMPLETED'},
                                    precursor_keys=[flow.root_task_key])
        for i in range(depth):
            child_task = flow.add_task(
                key=str(uuid4()),
                task={'status': 'COMPLETED', 'depth': depth},
                precursor_keys=[branch_root['key']]
            )
            branch_root = child_task
        flow.add_task(task={'status': 'PENDING'},
                      precursor_keys=[branch_root['key']])

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
        flow = Flow()
        flow.add_task(task={'status': 'COMPLETED'}, as_root=True)
        for i in range(3):
            flow.add_task(key=i, task={'status': 'RUNNING'},
                          precursor_keys=[flow.root_task_key])
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
        flow = Flow()
        flow.add_task(task={'status': 'COMPLETED'}, as_root=True)
        self.assertEqual(self.mocks['engine']['complete_flow'].call_count, 0)
        self.engine.tick_flow(flow)
        self.assertEqual(self.mocks['engine']['complete_flow'].call_count, 1)

class CompleteFlowTestCase(BaseTestCase):
    def test_sets_status_to_completed(self):
        flow = Flow()
        self.assertTrue(flow.status != 'COMPLETED')
        self.engine.complete_flow(flow=flow)
        self.assertTrue(flow.status == 'COMPLETED')

    def test_sets_output_for_single_tail_task(self):
        flow = Flow()
        tail_task = flow.add_task(task={'status': 'COMPLETED'},
                                  as_root=True)
        tail_task['output'] = 'some output'
        self.engine.complete_flow(flow=flow)
        self.assertEqual(flow.output, tail_task['output'])

    def test_sets_output_for_multiple_tail_tasks(self):
        flow = Flow()
        root_task = flow.add_task(task={'status': 'COMPLETED'},
                                  as_root=True)
        tail_tasks = []
        for i in range(3):
            tail_task = flow.add_task(task={'status': 'COMPLETED'},
                                      precursor_keys=[root_task['key']])
            tail_task['output'] = 'output_%s' % i
            tail_tasks.append(tail_task)
        self.engine.complete_flow(flow=flow)
        expected_output = {tail_task['key']: tail_task['output']
                           for tail_task in tail_tasks}
        self.assertEqual(flow.output, expected_output)

class StartTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.task = self.flow.add_task(task={})

    def test_processes_precomplete_actions(self):
        self.fail()

    def test_sets_status_to_running(self):
        self.engine.start_task(flow=self.flow, task=self.task)
        self.assertEqual(self.task['status'], 'RUNNING')

class CompleteTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()
        self.task = self.flow.add_task(task={})

    def test_processes_post_complete_actions(self):
        self.fail()


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
