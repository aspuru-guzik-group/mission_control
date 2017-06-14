from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock, patch

from .. import context_value_setter


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.setter = context_value_setter.ContextValueSetter()
        self.context = MagicMock()
        self.value_specs = [MagicMock() for i in range(3)]

    def setup_mocks(self, attrs=None):
        patchers = {attr: patch.object(context_value_setter, attr)
                    for attr in attrs}
        mocks = {}
        for key, patcher in patchers.items():
            self.addCleanup(patcher.stop)
            mocks[key] = patcher.start()
        return mocks

    def setup_setter_mocks(self, attrs=None):
        for attr in attrs: setattr(self.setter, attr, MagicMock())

class SeContextValuesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_setter_mocks(attrs=['set_context_value'])

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
        self.setup_setter_mocks(attrs=['get_value_for_value_spec'])

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

    def test_handles_from_source(self):
        self.value_spec['source'] = MagicMock()
        self.setter.get_value_from_dot_spec = MagicMock()
        result = self._get_value_for_value_spec()
        self.assertEqual(result,
                         self.setter.get_value_from_dot_spec.return_value)
        self.assertEqual(self.setter.get_value_from_dot_spec.call_args,
                         call(obj={'ctx': self.context},
                              dot_spec=self.value_spec['source']))

class TransformValueForValueSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.value = MagicMock()
        self.value_spec = defaultdict(MagicMock)
        self.context = MagicMock()

    def _transform(self, value=None, value_spec_updates=None):
        self.value_spec.update((value_spec_updates or {}))
        return self.setter.transform_value_for_value_spec(
            value=(value or self.value), value_spec=self.value_spec,
            context=self.context)

    def test_returns_value_if_no_transform_given(self):
        result = self._transform()
        self.assertEqual(result, self.value)

    @patch.object(context_value_setter, 'json')
    def test_json_dumps(self, _json):
        result = self._transform(value_spec_updates={'transform': 'json.dumps'})
        self.assertEqual(_json.dumps.call_args, call(self.value))
        self.assertEqual(result, _json.dumps.return_value)

    @patch.object(context_value_setter, 'json')
    def test_json_loads(self, _json):
        result = self._transform(value_spec_updates={'transform': 'json.loads'})
        self.assertEqual(_json.loads.call_args, call(self.value))
        self.assertEqual(result, _json.loads.return_value)

    def test_mapping_transform(self):
        self.setter.execute_mapping_transform = MagicMock()
        transform = {'type': 'mapping'}
        result = self._transform(value_spec_updates={'transform': transform})
        self.assertEqual(
            self.setter.execute_mapping_transform.call_args,
            call(value=self.value, transform=transform, context=self.context)
        )
        self.assertEqual(result,
                         self.setter.execute_mapping_transform.return_value)

class ExecuteMappingTransformTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.value = [MagicMock() for i in range(3)]
        self.transform = MagicMock()
        self.context = MagicMock()
        self.setup_setter_mocks(attrs=['map_item'])
        self.result = self.setter.execute_mapping_transform(
            value=self.value, transform=self.transform, context=self.context)

    def test_returns_mapped_items(self):
        self.assertEqual(
            self.setter.map_item.call_args_list,
            [call(item=item, idx=idx, items=self.value,
                  mapping_params=self.transform['params'], context=self.context)
             for idx, item in enumerate(self.value)]
        )
        self.assertEqual(
            self.result,
            [self.setter.map_item.return_value for item in self.value]
        )

class MapItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.item = MagicMock()
        self.items = MagicMock()
        self.idx = MagicMock()
        self.mapping_params = defaultdict(MagicMock, **{'skeleton': MagicMock(),
                                                        'wirings': MagicMock()})
        self.context = defaultdict(MagicMock)
        self.setup_setter_mocks(attrs=['set_context_values'])
        self.mocks = self.setup_mocks(attrs=['copy'])
        self.expected_copy = self.mocks['copy'].deepcopy.return_value
        self.expected_context = {
            **self.context, 'item': self.item, 'items': self.items,
            'idx': self.idx, 'skeleton': self.expected_copy
        }
        self.result = self.setter.map_item(item=self.item, idx=self.idx,
                                           items=self.items,
                                           mapping_params=self.mapping_params,
                                           context=self.context)

    def test_makes_deep_copy_of_the_skeleton(self):
        self.assertEqual(self.mocks['copy'].deepcopy.call_args,
                         call(self.mapping_params.get('skeleton')))

    def test_applies_wirings_with_expected_context(self):
        self.assertEqual(
            self.setter.set_context_values.call_args,
            call(value_specs=self.mapping_params.get('wirings'),
                 context=self.expected_context)
        )

    def test_returns_skeleton_copy(self):
        self.assertEqual(self.result, self.expected_copy)

class ExecuteMappingTransformE2ETestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.context = {'ctx_prop_{i}'.format(i=i): 'ctx_{i}_val'.format(i=i)
                        for i in range(3)}
        self.items = [
            {'prop_{i}'.format(i=i): 'item_{idx}_val_{i}'.format(idx=idx, i=i)
             for i in range(3)}
            for idx in range(3)
        ]
        self.skeleton = {'prop_{i}'.format(i=i): 'skel_{i}_val'.format(i=i)
                         for i in range(4)}
        self.transform = {
            'type': 'mapping',
            'params': {
                'skeleton': self.skeleton,
                'wirings': [
                    {'source': 'ctx.idx', 'dest': 'ctx.skeleton.idx'},
                    {'source': 'ctx.item.prop_0', 'dest': 'ctx.skeleton.prop_0'},
                    {'value': 'new_prop_1_val', 'dest': 'ctx.skeleton.prop_1'},
                    {'source': 'ctx.ctx_prop_2', 'dest': 'ctx.skeleton.prop_2'},
                ]
            }
        }
        self.result = self.setter.execute_mapping_transform(
            value=self.items, transform=self.transform, context=self.context)

    def test_produces_expected_items(self):
        expected_result = [
            {
                'idx': idx,
                'prop_0': item['prop_0'],
                'prop_1': 'new_prop_1_val',
                'prop_2': 'ctx_2_val',
                'prop_3': 'skel_3_val',
            }
            for idx, item in enumerate(self.items)
        ]
        self.assertEqual(self.result, expected_result)


