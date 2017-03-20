from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from ..set_values_task_handler import SetValuesTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task = defaultdict(MagicMock)
        self.job = defaultdict(MagicMock)
        self.task_handler = SetValuesTaskHandler()

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.set_values = MagicMock()

    def _do_initial_tick(self):
        self.task_handler.initial_tick(task=self.task, job=self.job)

    def test_calls_set_values(self):
        self._do_initial_tick()
        self.assertEqual(
            self.task_handler.set_values.call_args,
            call(value_specs=self.task['params']['value_specs'],
                 job=self.job))

    def test_sets_task_status_to_completed(self):
        self._do_initial_tick()
        self.assertEqual(self.task['status'], 'COMPLETED')

class SetValuesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.set_context_value = MagicMock()
        self.task_handler.get_context = MagicMock()
        self.expected_context = self.task_handler.get_context.return_value
        self.value_specs = [MagicMock() for i in range(3)]

    def test_calls_set_context_value_for_value_specs(self):
        self.task_handler.set_values(value_specs=self.value_specs, job=self.job)
        self.assertEqual(
            self.task_handler.set_context_value.call_args_list,
            [call(value_spec=value_spec, context=self.expected_context)
             for value_spec in self.value_specs]
        )

class GetContext(BaseTestCase):
    def test_context_has_keyed_tasks(self):
        self.job['tasks'] = [MagicMock() for i in range(3)]
        context = self.task_handler.get_context(job=self.job)
        expected_keyed_tasks = {task.get('key'): task
                                for task in self.job['tasks']}
        self.assertEqual(context['keyed_tasks'], expected_keyed_tasks)

    def test_context_has_job(self):
        context = self.task_handler.get_context(job=self.job)
        self.assertEqual(context['job'], self.job)

class SetContextValueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.target_dot_spec = 'foo.bar'
        self.context = {}
        self.src_key = 'some_key'
        self.context[self.src_key] = 'src_value'
        self.value_spec = {
            'target': self.target_dot_spec,
            'value': '{{%s}}' % self.src_key
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
