from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from .. import task_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task_handlers = defaultdict(MagicMock)
        self.task = defaultdict(MagicMock)
        self.job = defaultdict(MagicMock)
        self.task_runner = task_runner.TaskRunner(
            task_handlers=self.task_handlers)

class TickTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_runner.get_task_handler = MagicMock()

    def test_dispatches_to_task_handler(self):
        self.task_runner.tick_task(task=self.task, job=self.job)
        expected_task_handler = self.task_runner.get_task_handler.return_value
        self.assertEqual(expected_task_handler.tick_task.call_args,
                         call(task=self.task, job=self.job))

class GetTaskHandlerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task['type'] = 'some_type'
        self.task_runner.task_handlers[self.task['type']] = MagicMock()
    
    def test_gets_handler(self):
        result = self.task_runner.get_task_handler(task=self.task, job=self.job)
        self.assertEqual(result,
                         self.task_runner.task_handlers[self.task['type']])
