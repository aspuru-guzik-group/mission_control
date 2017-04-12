from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import set_values_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task = defaultdict(MagicMock)
        self.task_context = defaultdict(MagicMock)
        self.task_handler = set_values_task_handler.SetValuesTaskHandler()

@patch.object(set_values_task_handler, 'set_context_values')
class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.set_values = MagicMock()

    def _do_initial_tick(self):
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_calls_set_context_values(self, mock_set_context_values):
        self._do_initial_tick()
        self.assertEqual(
            mock_set_context_values.call_args,
            call(value_specs=self.task['task_params']['value_specs'],
                 context=self.task_context))

    def test_sets_task_status_to_completed(self, *args):
        self._do_initial_tick()
        self.assertEqual(self.task['status'], 'COMPLETED')
