import unittest
from unittest.mock import call, MagicMock

from .. import context_value_setter


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.setter = context_value_setter.ContextValueSetter()
        self.context = MagicMock()
        self.value_specs = [MagicMock() for i in range(3)]

class SeContextValuesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setter.set_context_value = MagicMock()

    def test_calls_set_context_value_for_value_specs(self):
        self.setter.set_context_values(value_specs=self.value_specs,
                                       context=self.context)
        self.assertEqual(
            self.setter.set_context_value.call_args_list,
            [call(value_spec=value_spec, context=self.context)
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
        self.setter.set_context_value(value_spec=self.value_spec,
                                      context=self.context)
        self.assertEqual(self.context['foo']['bar'],
                         self.context.get(self.src_key))

    def test_sets_raw_value_on_context(self):
        self.value_spec['is_raw_value'] = True
        self.setter.set_context_value(value_spec=self.value_spec,
                                      context=self.context)
        self.assertEqual(self.context['foo']['bar'], self.value_spec['value'])

class GetValueForValueSpecTestCase(BaseTestCase):
    def test_handles_from_value(self):
        self.fail()

    def test_handles_template_value(self):
        self.fail()

    def test_handles_context_value(self):
        self.fail()

    def test_handles_json_value(self):
        self.fail()
