from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from ..execute_job_task_handler import ExecuteJobTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task = self.generate_task()
        self.job = defaultdict(MagicMock)
        self.task_context = defaultdict(MagicMock)
        self.task_context['job'] = self.job
        self.execution_client = MagicMock()
        self.task_handler = ExecuteJobTaskHandler(
            execution_client=self.execution_client
        )

    def generate_task(self):
        task = defaultdict(MagicMock)
        task['data'] = defaultdict(MagicMock)
        return task

class InitialTickTestCase(BaseTestCase):
    def _do_initial_tick(self):
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_starts_execution(self):
        self._do_initial_tick()
        self.assertEqual(
            self.task_handler.execution_client.start_execution.call_args,
            call(submission=self.task['task_params']['submission']))

    def test_stores_execution_meta(self):
        self._do_initial_tick()
        self.assertEqual(
            self.task['data']['execution_meta'],
            self.task_handler.execution_client.start_execution.return_value
        )

class IntermediateTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.handle_execution_state = MagicMock()

    def _do_intermediate_tick(self):
        self.task_handler.intermediate_tick(task=self.task,
                                            task_context=self.task_context)

    def test_calls_get_execution_state_with_execution_meta(self):
        self._do_intermediate_tick()
        self.assertEqual(
            self.task_handler.execution_client.get_execution_state.call_args,
            call(execution_meta=self.task['data']['execution_meta']))

    def test_updates_job_execution_state(self):
        self._do_intermediate_tick()
        self.assertEqual(
            self.task['data']['execution_state'],
            self.task_handler.execution_client.get_execution_state.return_value)

    def test_calls_handle_execution_state(self):
        self._do_intermediate_tick()
        self.assertEqual(
            self.task_handler.handle_execution_state.call_args,
            call(execution_state=self.task_handler.execution_client\
                 .get_execution_state.return_value,
                 task=self.task,
                 job=self.task_context['job'])
        )

class HandleExecutionStateTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.handle_completed_execution_state = MagicMock()
        self.task_handler.handle_failed_execution_state = MagicMock()

    def _do_handle_execution_state(self, execution_state=None):
        self.task_handler.handle_execution_state(
            execution_state=execution_state, task=self.task, job=self.job)

    def test_handles_completed_execution(self):
        completed_execution_state = {'run_status': 'COMPLETED'}
        self._do_handle_execution_state(
            execution_state=completed_execution_state)
        self.assertEqual(
            self.task_handler.handle_completed_execution_state.call_args,
            call(execution_state=completed_execution_state,
                 task=self.task,
                 job=self.job)
        )

    def test_handles_failed_execution(self):
        failed_execution_state = {'run_status': 'FAILED'}
        self._do_handle_execution_state(execution_state=failed_execution_state)
        self.assertEqual(
            self.task_handler.handle_failed_execution_state.call_args,
            call(execution_state=failed_execution_state,
                 task=self.task,
                 job=self.job)
        )

class HandleCompletedExecutionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.completed_execution_state = MagicMock()

    def _do_handle_completed_execution_state(self):
        self.task_handler.handle_completed_execution_state(
            execution_state=self.completed_execution_state,
            task=self.task,
            job=self.job)

    def test_saves_artifact_to_task_data(self):
        self._do_handle_completed_execution_state()
        self.assertEqual(self.task['data']['artifact'],
                         self.completed_execution_state.get('artifact'))

    def test_saves_stdout_to_task_data(self):
        self._do_handle_completed_execution_state()
        self.assertEqual(self.task['data']['stdout'],
                         self.completed_execution_state.get('stdout'))

    def test_sets_status(self):
        self._do_handle_completed_execution_state()
        self.assertEqual(self.task['status'], 'COMPLETED')

class HandleFailedExecutionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.failed_execution_state = {'error': 'some_error'}

    def _do_handle_failed_execution_state(self):
        self.task_handler.handle_failed_execution_state(
            execution_state=self.failed_execution_state,
            task=self.task,
            job=self.job)

    def test_saves_failure_to_task(self):
        self._do_handle_failed_execution_state()
        self.assertEqual(self.task['data']['error'],
                         self.failed_execution_state['error'])

    def test_sets_status(self):
        self._do_handle_failed_execution_state()
        self.assertEqual(self.task['status'], 'FAILED')
