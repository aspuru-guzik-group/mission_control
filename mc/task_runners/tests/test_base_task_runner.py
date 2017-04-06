from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from .. import base_task_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.tasks = self.generate_tasks()
        self.get_tasks = lambda: self.tasks
        self.task_context = MagicMock()
        self.get_task_context = lambda: self.task_context
        self.task_handler = MagicMock()
        self.task_runner = base_task_runner.BaseTaskRunner(
            get_tasks=self.get_tasks,
            task_handler=self.task_handler,
            get_task_context=self.get_task_context
        )

    def generate_tasks(self):
        tasks = [defaultdict(MagicMock) for i in range(3)]
        return tasks

class TickTasksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.completing_tasks = [{'key': i, 'finish': 0} for i in range(3)]
        self.noncompleting_tasks = [{'key': i, 'finish': 1} for i in range(3)]
        def mock_tick_task(task=None):
            if task in self.completing_tasks: task['status'] = 'COMPLETED'
            else: task['status'] = 'RUNNING'
        self.task_runner.tick_task = MagicMock(side_effect=mock_tick_task)

    def test_continues_until_first_noncompleting_task(self):
        self.tasks = self.completing_tasks + self.noncompleting_tasks
        self.task_runner.tick_tasks()
        expected_tick_task_call_args_list = (
            [call(task=task) for task in self.completing_tasks] +
            [call(task=self.noncompleting_tasks[0])]
        )
        self.assertEqual(self.task_runner.tick_task.call_args_list,
                         expected_tick_task_call_args_list)

    def test_returns_running_status_if_has_incomplete_tasks(self):
        self.tasks = self.noncompleting_tasks
        tasks_status = self.task_runner.tick_tasks()
        self.assertEqual(tasks_status, 'RUNNING')

    def test_returns_completed_status_if_no_incomplete_tasks(self):
        self.tasks = self.completing_tasks
        tasks_status = self.task_runner.tick_tasks()
        self.assertEqual(tasks_status, 'COMPLETED')
        self.task_runner.tasks = []
        tasks_status = self.task_runner.tick_tasks()
        self.assertEqual(tasks_status, 'COMPLETED')

class TickTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_runner.get_decorated_task_context = MagicMock()
        self.expected_task_context = \
                self.task_runner.get_decorated_task_context.return_value

    def test_dispatches_to_task_handler(self):
        task = self.tasks[0]
        self.task_runner.tick_task(task=task)
        self.assertEqual(
            self.task_runner.task_handler.tick_task.call_args,
            call(task=task, task_context=self.expected_task_context)
        )

class DecorateContextTestCase(BaseTestCase):
    def test_decorated_task_context_has_keyed_tasks(self):
        self.task_context = {'key_%s' % i : 'value_%s' % i for i in range(3)}
        expected_keyed_tasks = {task['task_key']: task
                                for task in self.task_runner.get_tasks()}
        expected_decorated_task_context = {
            'tasks': expected_keyed_tasks,
            **self.task_context,
        }
        self.assertEqual(self.task_runner.get_decorated_task_context(),
                         expected_decorated_task_context)
