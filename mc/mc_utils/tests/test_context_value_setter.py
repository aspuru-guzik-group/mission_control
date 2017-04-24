from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock, patch

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
        self.value_spec = {'dest': self.dest_dot_spec}
        self.setter.get_value_for_value_spec = MagicMock()

    def test_sets_dest_value_on_context(self):
        self.setter.set_context_value(value_spec=self.value_spec,
                                      context=self.context)
        self.assertEqual(self.context['foo']['bar'],
                         self.setter.get_value_for_value_spec.return_value)

class GetValueForValueSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.value_spec = defaultdict(MagicMock, **{'dest': MagicMock()})
        self.context = MagicMock()

    def _get_value_for_value_spec(self):
        return self.setter.get_value_for_value_spec(value_spec=self.value_spec,
                                                    context=self.context)

    def test_handles_from_value(self):
        self.value_spec['value'] = MagicMock()
        result = self._get_value_for_value_spec()
        self.assertEqual(result, self.value_spec['value'])

    def test_handles_template_value(self):
        self.setter.render_template = MagicMock()
        self.value_spec['template'] = MagicMock()
        result = self._get_value_for_value_spec()
        self.assertEqual(result, self.setter.render_template.return_value)
        self.assertEqual(self.setter.render_template.call_args,
                         call(template=self.value_spec['template'],
                              context=self.context))

    @patch.object(context_value_setter, 'json')
    def test_handles_json_value(self, mock_json):
        self.value_spec['from_json'] = True
        self.setter.get_value_from_dot_spec = MagicMock()
        result = self._get_value_for_value_spec()
        self.assertEqual(result, mock_json.loads.return_value)
        self.assertEqual(mock_json.loads.call_args,
                         call(self.setter.get_value_from_dot_spec.return_value))
        self.assertEqual(self.setter.get_value_from_dot_spec.call_args,
                         call(obj={'ctx': self.context},
                              dot_spec=self.value_spec['source']))

    def test_handles_dot_spec(self):
        self.setter.get_value_from_dot_spec = MagicMock()
        result = self._get_value_for_value_spec()
        self.assertEqual(result,
                         self.setter.get_value_from_dot_spec.return_value)
        self.assertEqual(self.setter.get_value_from_dot_spec.call_args,
                         call(obj={'ctx': self.context},
                              dot_spec=self.value_spec['source']))
