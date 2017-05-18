from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import wire_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.wirings = [MagicMock() for i in range(3)]
        self.task = defaultdict(MagicMock, )
        self.task_context = defaultdict(MagicMock, **{
            'task_params': {'wirings': self.wirings}
        })
        self.task_handler = wire_task_handler.WireTaskHandler()

@patch.object(wire_task_handler, 'context_value_setter')
class InitialTickTestCase(BaseTestCase):
    def _do_initial_tick(self):
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_calls_set_context_value(self, _context_value_setter):
        self._do_initial_tick()
        expected_set_value_call_args_list = [
            call(value_spec=wiring, task_context=self.task_context)
            for wiring in self.task['task_params']['wirings']
        ]
        self.assertEqual(_context_value_setter.set_value.call_args_list,
                         expected_set_value_call_args_list)

    def test_sets_task_status_to_completed(self, *args):
        self._do_initial_tick()
        self.assertEqual(self.task['status'], 'COMPLETED')
