import collections
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

class FlowSpecToFlowTestCase(BaseTestCase):
    @patch.object(flow_engine, 'Flow')
    def test_dispatches_to_flow(self, _Flow):
        flow_spec = MagicMock()
        result = flow_engine.FlowEngine.flow_spec_to_flow(flow_spec=flow_spec)
        self.assertEqual(_Flow.from_flow_spec.call_args,
                         call(flow_spec=flow_spec))
        self.assertEqual(result, _Flow.from_flow_spec.return_value)

class FlowDictToFlowTestCase(BaseTestCase):
    @patch.object(flow_engine, 'Flow')
    def test_dispatches_to_Flow(self, _Flow):
        flow_dict = MagicMock()
        result =flow_engine.FlowEngine.flow_dict_to_flow(flow_dict=flow_dict)
        self.assertEqual(_Flow.from_flow_dict.call_args,
                         call(flow_dict=flow_dict))
        self.assertEqual(result, _Flow.from_flow_dict.return_value)

class FlowToFlowDictTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = MagicMock()
        self.result = self.engine.flow_to_flow_dict(flow=self.flow)

    def test_dispatches_to_flow(self):
        self.assertEqual(self.result, self.flow.to_flow_dict.return_value)

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

    def test_raises_if_flow_fails(self):
        self._setup_mock_tick(flow_end_status='FAILED')
        with self.assertRaises(self.engine.FlowError):
            self.engine.run_flow(flow=self.flow)
            self.assertEqual(self.flow.status, 'FAILED')

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {'engine': patch.multiple(
            self.engine, start_task=DEFAULT, complete_flow=DEFAULT)}
        self.mocks = self.start_patchers(patchers=self.patchers)
        def mock_tick_task(*args, flow=None, task=None, task_ctx=None):
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
        self.task = collections.defaultdict(MagicMock)
        self.task_ctx = collections.defaultdict(MagicMock)
        self.engine.complete_task = MagicMock()
        self.engine.tick_proxying_task = MagicMock()

    def _tick_task(self):
        self.engine.tick_task(flow=self.flow, task=self.task,
                              task_ctx=self.task_ctx)

    def test_calls_tick_proxying_task_if_task_is_proxy_task(self):
        self.task['proxied_task'] = MagicMock()
        self._tick_task()
        self.assertEqual(self.engine.tick_proxying_task.call_args,
                         call(proxying_task=self.task, flow=self.flow,
                              task_ctx=self.task_ctx))

    def test_calls_task_handler(self):
        self._tick_task()
        expected_task_ctx = {**self.task_ctx, 'task': self.task,
                             'flow': self.flow}
        self.assertEqual(self.engine.task_handler.tick_task.call_args,
                         call(task_ctx=expected_task_ctx))

    def test_calls_complete_task_if_task_completes(self):
        def mock_tick_task(*args, **kwargs): self.task['status'] = 'COMPLETED'
        self.engine.task_handler.tick_task.side_effect = mock_tick_task
        self._tick_task()
        self.assertEqual(self.engine.complete_task.call_args,
                         call(flow=self.flow, task=self.task))

class TickProxyingTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.proxied_task = MagicMock()
        self.proxying_task = collections.defaultdict(MagicMock, {
            'proxied_task': self.proxied_task
        })
        self.flow = MagicMock()
        self.task_ctx = MagicMock()
        self.engine.tick_task = MagicMock()
        self.engine.tick_proxying_task(
            proxying_task=self.proxying_task, flow=self.flow,
            task_ctx=self.task_ctx)

    def test_calls_tick_task_on_proxied_task(self):
        self.assertEqual(self.engine.tick_task.call_args,
                         call(task=self.proxied_task, flow=self.flow,
                              task_ctx=self.task_ctx))

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
