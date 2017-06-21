import types
import unittest
from unittest.mock import call, MagicMock, patch

from .. import dot_spec_loader

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.loader = dot_spec_loader.DotSpecLoader()

class _LoadFromDotSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        importlib_patcher = patch.object(dot_spec_loader, 'importlib')
        self.addCleanup(importlib_patcher.stop)
        self.mock_importlib = importlib_patcher.start()
        self.loader._get_obj_value_from_dot_spec = MagicMock()
        self.module_part = '1.2.3'
        self.obj_part = 'a.b.c'
        self.dot_spec = self.module_part + ':' + self.obj_part
        self.result = self.loader._load_from_dot_spec(dot_spec=self.dot_spec)

    def test_imports_module(self):
        self.assertEqual(self.mock_importlib.import_module.call_args,
                         call(self.module_part))

    def test_gets_obj_value_from_module(self):
        self.assertEqual(
            self.loader._get_obj_value_from_dot_spec.call_args,
            call(obj=self.mock_importlib.import_module.return_value,
                 dot_spec=self.obj_part)
        )
        self.assertEqual(self.result,
                         self.loader._get_obj_value_from_dot_spec.return_value)

class _GetObjValueFromDotSpecTestCase(BaseTestCase):
    def test_calls_get_attr_or_item(self):
        self.loader._get_attr_or_item = MagicMock()
        obj = MagicMock()
        path_elements = [str(i) for i in range(5)]
        dot_spec = '.'.join(path_elements)
        result = self.loader._get_obj_value_from_dot_spec(obj, dot_spec)
        expected_call_args_list = []
        cursor = obj
        for path_element in path_elements:
            expected_call_args_list.append(call(obj=cursor, key=path_element))
            next_cursor = self.loader._get_attr_or_item.return_value
            cursor = next_cursor
        self.assertEqual(self.loader._get_attr_or_item.call_args_list,
                         expected_call_args_list)
        self.assertEqual(result, self.loader._get_attr_or_item.return_value)

class _GetObjValueFromDotSpecTestCaseE2E(BaseTestCase):
    def test_handles_dict_views(self):
        obj = {i: MagicMock() for i in range(3)}
        result = self.loader._get_obj_value_from_dot_spec(
            obj=obj, dot_spec='values().0')
        self.assertEqual(result, list(obj.values())[0])

class _GetAttrOrItemTestCase(BaseTestCase):
    def _get(self, obj, key):
        return self.loader._get_attr_or_item(obj=obj, key=key)

    def test_gets_from_sequence(self):
        seq = [MagicMock() for i in range(3)]
        key = 0
        self.assertEqual(self._get(seq, key),  seq[key])

    def test_gets_from_mapping(self):
        mapping = {'key_%s' % i: MagicMock() for i in range(3)}
        key = list(mapping.keys())[1]
        self.assertEqual(self._get(mapping, key),  mapping[key])

    def test_gets_from_callable(self):
        mapping = {'key_%s' % i: MagicMock() for i in range(3)}
        self.assertEqual(self._get(mapping, 'items()'),  mapping.items())

    def test_gets_from_attr(self):
        obj = types.SimpleNamespace()
        attrs = ['attr_%s' % i for i in range(3)]
        for attr in attrs: setattr(obj, attr, MagicMock())
        key = attrs[1]
        self.assertEqual(self._get(obj, key),  getattr(obj, key))
