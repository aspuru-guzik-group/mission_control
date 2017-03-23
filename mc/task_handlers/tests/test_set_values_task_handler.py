from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from ..set_values_task_handler import SetValuesTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task = defaultdict(MagicMock)
        self.task_context = defaultdict(MagicMock)
        self.task_handler = SetValuesTaskHandler()

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.set_values = MagicMock()

    def _do_initial_tick(self):
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_calls_set_values(self):
        self._do_initial_tick()
        self.assertEqual(
            self.task_handler.set_values.call_args,
            call(value_specs=self.task['params']['value_specs'],
                 task_context=self.task_context))

    def test_sets_task_status_to_completed(self):
        self._do_initial_tick()
        self.assertEqual(self.task['status'], 'COMPLETED')

class SetValuesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.set_context_value = MagicMock()
        self.value_specs = [MagicMock() for i in range(3)]

    def test_calls_set_context_value_for_value_specs(self):
        self.task_handler.set_values(value_specs=self.value_specs,
                                     task_context=self.task_context)
        self.assertEqual(
            self.task_handler.set_context_value.call_args_list,
            [call(value_spec=value_spec, context=self.task_context)
             for value_spec in self.value_specs]
        )

class SetContextValueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dest_dot_spec = 'ctx.foo.bar'
        self.context = {}
        self.src_key = 'some_key'
        self.context[self.src_key] = 'src_value'
        self.value_spec = {
            'dest': self.dest_dot_spec,
            'value': '{{ctx.%s}}' % self.src_key
        }

    def test_sets_transformed_value_on_context(self):
        self.task_handler.set_context_value(
            value_spec=self.value_spec,
            context=self.context
        )
        self.assertEqual(self.context['foo']['bar'],
                         self.context.get(self.src_key))

    def test_sets_raw_value_on_context(self):
        self.value_spec['is_raw_value'] = True
        self.task_handler.set_context_value(
            value_spec=self.value_spec,
            context=self.context
        )
        self.assertEqual(self.context['foo']['bar'], self.value_spec['value'])
