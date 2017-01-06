import unittest
from unittest.mock import call, DEFAULT, patch, Mock
from uuid import uuid4

from ..flow_engine import FlowEngine
from ..flow import Flow, BaseTask


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
        self.setup_task_classes()

    def setup_task_classes(self):
        class BaseTask(object):
            def __init__(self, *args, key=None, status=None, ctx=None, **kwargs):
                self.key = key
                self.status = status
                self.ctx = ctx

        self.task_classes = {name: type(name, (BaseTask,), {})
                             for name in ['Type1', 'Type2', 'Type3']}
        for name, task_class in self.task_classes.items():
            self.engine.register_task_class(task_class=task_class)

        self.serialized_flow = {
            'data': 'some data',
            'input': 'some input',
            'output': 'some output',
            'tasks': [
                {'key': 'a', 'task_type': 'Type1', 'status': 'COMPLETED'},
                {'key': 'b', 'task_type': 'Type2', 'status': 'COMPLETED'},
                {'key': 'c', 'task_type': 'Type1', 'status': 'RUNNING'},
                {'key': 'd', 'task_type' : 'Type2', 'status': 'RUNNING'},
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
        expected_tasks = {}
        for serialized_task in self.serialized_flow['tasks']:
            expected_task = self.engine.deserialize_task(serialized_task)
            expected_tasks[expected_task.key] = expected_task

        def summarize_tasks(tasks):
            return {
                task.key: {
                    'key': task.key,
                    'task_type': type(task),
                    'status': task.status,
                    'ctx': task.ctx,
                }
                for task in tasks.values()
            }

        self.assertEqual(summarize_tasks(self.flow.tasks),
                         summarize_tasks(expected_tasks))

    def test_has_expected_root_task(self):
        self.assertEqual(self.flow.root_task, self.flow.tasks['a'])

    def test_has_expected_edges(self):
        expected_edges = {}
        for serialized_edge in self.serialized_flow['edges']:
            src_key, dest_key = (serialized_edge['src_key'],
                               serialized_edge['dest_key'])
            expected_edges[(src_key, dest_key)] = {
                'src': self.flow.tasks[src_key],
                'dest': self.flow.tasks[dest_key]
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

        flow.input = 'some input'
        flow.output = 'some output'
        flow.data = 'some data'

        class BaseTask(object):
            task_type = 'BaseTask'
            def __init__(self, key=None, status=None, state=None, **kwargs):
                self.key = key
                self.status = status
                self.state = state

        tasks = [BaseTask(key=i, status='s_%s' % i, state=i) for i in range(3)]
        flow.add_tasks(tasks=tasks)
        flow.root_task = tasks[0]
        edges = [{'src': tasks[0], 'dest': tasks[1]},
                 {'src': tasks[1], 'dest': tasks[2]}]
        flow.add_edges(edges=edges)
        flow.state = 'flow_state'
        flow.status = 'flow_status'
        return flow

    def test_serializes_flow(self):
        serialization = self.engine.serialize_flow(self.flow)
        expected_serialization = {
            'data': self.flow.data,
            'input': self.flow.input,
            'output': self.flow.output,
            'tasks': [self.engine.serialize_task(task=task)
                      for task in self.flow.tasks.values()],
            'root_task_key': self.flow.root_task.key,
            'edges': [{'src_key': edge['src'].key, 'dest_key': edge['dest'].key}
                      for edge in self.flow.edges.values()],
            'state': self.flow.state,
            'status': self.flow.status,
        }
        self.assertEqual(serialization, expected_serialization)

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {'engine': patch.multiple(
            self.engine, start_task=DEFAULT, complete_flow=DEFAULT)}
        self.mocks = self.start_patchers(patchers=self.patchers)

    def test_starts_nearest_pending_tasks(self):
        self.maxDiff = None
        flow = self.generate_flow_with_pending_successors()
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

    def generate_flow_with_pending_successors(self):
        flow = Flow()
        root_task = Mock()
        root_task.configure_mock(key='ROOT', status='COMPLETED')
        flow.add_task(root_task, as_root=True)
        self.add_branch_with_pending_leaf_to_flow(flow, depth=3)
        self.add_branch_with_pending_leaf_to_flow(flow, depth=2)
        return flow

    def add_branch_with_pending_leaf_to_flow(self, flow, depth=3):
        branch_root = Mock()
        branch_root.status = 'COMPLETED'
        flow.add_task(task=branch_root, precursor=flow.root_task)
        for i in range(depth):
            child_task = Mock()
            child_task.key = str(uuid4())
            child_task.status = 'COMPLETED'
            child_task.depth = depth
            flow.add_task(task=child_task, precursor=branch_root)
            branch_root = child_task
        branch_leaf = Mock()
        branch_leaf.status = 'PENDING'
        flow.add_task(task=branch_leaf, precursor=branch_root)

    def test_ticks_running_tasks(self):
        flow = self.generate_flow_with_running_tasks()
        running_tasks = flow.get_tasks_by_status(status='RUNNING')
        expected_task_summaries_before_tick = {
            task.key: self.summarize_task(task, overrides={'tick_count': 0,
                                                          'status': 'RUNNING'})
            for task in running_tasks
        }
        self.assertEqual(self.summarize_tasks(running_tasks),
                         expected_task_summaries_before_tick)
        self.engine.tick_flow(flow)
        expected_task_summaries_after_tick = {
            task.key: self.summarize_task(task, overrides={'tick_count': 1,
                                                          'status': 'RUNNING'})
            for task in running_tasks
        }
        self.assertEqual(self.summarize_tasks(running_tasks),
                         expected_task_summaries_after_tick)

    def generate_flow_with_running_tasks(self):
        flow = Flow()
        flow.add_task(BaseTask(status='COMPLETED'), as_root=True)
        for i in range(3):
            running_task = Mock()
            running_task.configure_mock(key=i, status='RUNNING')
            flow.add_task(running_task, precursor=flow.root_task)
        return flow

    def summarize_tasks(self, tasks):
        return {task.key: self.summarize_task(task) for task in tasks}

    def summarize_task(self, task, overrides=None):
        if not overrides: overrides = {}
        return {'key': task.key, 'tick_count': task.tick.call_count,
                'status': task.status, **overrides}

    def test_calls_complete_flow_if_no_incomplete_tasks(self):
        flow = Flow()
        flow.add_task(task=BaseTask(status='COMPLETED'), as_root=True)
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
        tail_task = flow.add_task(task=BaseTask(status='COMPLETED'),
                                  as_root=True)
        tail_task.output = 'some output'
        self.engine.complete_flow(flow=flow)
        self.assertEqual(flow.output, tail_task.output)

    def test_sets_output_for_multiple_tail_tasks(self):
        flow = Flow()
        root_task = flow.add_task(task=BaseTask(status='COMPLETED'),
                                  as_root=True)
        tail_tasks = []
        for i in range(3):
            tail_task = flow.add_task(task=BaseTask(status='COMPLETED'),
                                      precursor=root_task)
            tail_task.output = 'output_%s' % i
            tail_tasks.append(tail_task)
        self.engine.complete_flow(flow=flow)
        expected_output = {tail_task.key: tail_task.output
                          for tail_task in tail_tasks}
        self.assertEqual(flow.output, expected_output)

class StartTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = Flow()

    def test_sets_status_to_running(self):
        task = self.flow.add_task(task=Mock())
        self.engine.start_task(flow=self.flow, task=task)
        self.assertEqual(task.status, 'RUNNING')

    def test_sets_input_for_multiple_precursors(self):
        precursors = [self.flow.add_task(task=Mock()) for i in range(3)]
        successor = self.flow.add_task(task=Mock(), precursor=precursors)
        successor.static_input = None
        self.engine.start_task(flow=self.flow, task=successor)
        expected_input = [precursor.output for precursor in precursors]
        self.assertEqual(set(successor.input), set(expected_input))

    def test_sets_input_for_single_precursor(self):
        precursor = self.flow.add_task(task=Mock())
        successor = self.flow.add_task(task=Mock(), precursor=precursor)
        successor.static_input = None
        self.engine.start_task(flow=self.flow, task=successor)
        expected_input = precursor.output
        self.assertEqual(successor.input, expected_input)

    def test_uses_static_input_if_provided(self):
        precursor = self.flow.add_task(task=Mock())
        successor = self.flow.add_task(task=Mock(), precursor=precursor)
        successor.static_input = Mock()
        self.engine.start_task(flow=self.flow, task=successor)
        expected_input = successor.static_input
        self.assertEqual(successor.input, expected_input)

class GenerateFlowTestCase(BaseTestCase):
    def test_generates_flow_from_spec(self):
        MockFlowGenerator = Mock()
        self.engine.register_flow_generator(MockFlowGenerator,
                                            test_fn=lambda flow_spec: True)
        flow_spec = Mock()
        flow = self.engine.generate_flow(flow_spec=flow_spec)
        self.assertEqual(MockFlowGenerator.generate_flow.call_args,
                         call(flow_spec=flow_spec))
        self.assertEqual(flow, MockFlowGenerator.generate_flow.return_value)

if __name__ == '__main__':
    unittest.main()
