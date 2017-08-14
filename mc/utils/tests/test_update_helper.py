import unittest
from types import SimpleNamespace
from .. import update_helper


class BaseTestCase(unittest.TestCase):
    pass


class DictGetNestedTestCase(BaseTestCase):
    def test_gets_nested_handle(self):
        d = {'a': {'b': {'c': 1}}}
        handle = update_helper.update_helper._get_nested_handle(d, 'a.b.c')
        key, accessors = handle['key'], handle['accessors']
        self.assertEqual(accessors['getter'](key), 1)
        accessors['setter'](key, 2)
        self.assertEqual(d['a']['b']['c'], 2)
        accessors['deleter'](key)
        self.assertEqual(d['a']['b'], {})

    def test_creates_nonexistent_keys(self):
        d = {'a': {}}
        update_helper.update_helper._get_nested_handle(d, 'a.b.c')
        self.assertEqual(d, {'a': {'b': {}}})


class ObjGetNestedTestCase(BaseTestCase):
    def test_gets_nested(self):
        d = SimpleNamespace()
        parent = d
        for token in ['a', 'b', 'c']:
            setattr(parent, token, SimpleNamespace())
            parent = getattr(parent, token)
        d.a.b.c = 1
        handle = update_helper.update_helper._get_nested_handle(d, 'a.b.c')
        key, accessors = handle['key'], handle['accessors']
        self.assertEqual(accessors['getter'](key), 1)
        accessors['setter'](key, 2)
        self.assertEqual(d.a.b.c, 2)
        accessors['deleter'](key)
        self.assertFalse(hasattr(d.a.b, 'c'))

    def test_creates_nonexistent_keys(self):
        d = SimpleNamespace()
        update_helper.update_helper._get_nested_handle(d, 'a.b.c.d')
        self.assertEqual(d.a, {'b': {'c': {}}})


class AddTestCase(BaseTestCase):
    def test_dict(self):
        d = {'a': 1}
        update_helper.update(d, [('a', '$add', 3)])
        self.assertEqual(d['a'], 4)

    def test_obj(self):
        d = SimpleNamespace(a=1)
        update_helper.update(d, [('a', '$add', 3)])
        self.assertEqual(d.a, 4)


class MulTestCase(BaseTestCase):
    def test_multiplies(self):
        d = {'a': 1}
        update_helper.update(d, [('a', '$mul', 3)])
        self.assertEqual(d['a'], 3)


class SetTestCase(BaseTestCase):
    def test_sets(self):
        d = {'a': 1}
        update_helper.update(d, [('a', '$set', 2)])
        self.assertEqual(d['a'], 2)


class UnsetTestCase(BaseTestCase):
    def test_deletes(self):
        d = {'a': 1}
        update_helper.update(d, [('a', '$unset')])
        self.assertTrue('a' not in d)


class ShiftTestCase(BaseTestCase):
    def test_shifts(self):
        d = {'a': [1, 2, 3]}
        update_helper.update(d, [('a', '$shift')])
        self.assertEqual(d['a'], [2, 3])


class UnshiftTestCase(BaseTestCase):
    def test_unshifts(self):
        d = {'a': [1, 2, 3]}
        update_helper.update(d, [('a', '$unshift', [7, 8, 9])])
        self.assertEqual(d['a'], [7, 8, 9, 1, 2, 3])


class SpliceTestCase(BaseTestCase):
    def test_splices(self):
        d = {'a': [1, 2, 3]}
        update_helper.update(d, [('a', '$splice', {'start': 1,
                                                   'delete_count': 2,
                                                   'new_items': [7, 8, 9]})])
        self.assertEqual(d['a'], [1, 7, 8, 9])


class MergeTestCase(BaseTestCase):
    def test_merges(self):
        d = {'a': {'b': 'b.1', 'c': 'c.1'}}
        update_helper.update(d, [('a', '$merge', {'b': 'b.2', 'd': 'd.2'})])
        self.assertEqual(d['a'], {'b': 'b.2', 'c': 'c.1', 'd': 'd.2'})


class PushTestCase(BaseTestCase):
    def test_pushes(self):
        d = {'a': [1, 2, 3]}
        update_helper.update(d, [('a', '$push', [7, 8, 9])])
        self.assertEqual(d['a'], [1, 2, 3, 7, 8, 9])


class PopTestCase(BaseTestCase):
    def test_pops(self):
        d = {'a': [1, 2, 3]}
        update_helper.update(d, [('a', '$pop')])
        self.assertEqual(d['a'], [1, 2])


class OmitTestCase(BaseTestCase):
    def test_omits(self):
        d = {'a': {'b': 'b', 'c': 'c', 'd': 'd'}}
        update_helper.update(d, [('a', '$omit', ['b', 'd'])])
        self.assertEqual(d['a'], {'c': 'c'})


class AddToSetTestCase(BaseTestCase):
    def test_adds_to_set(self):
        d = {'a': [1, 2, 3]}
        update_helper.update(d, [('a', '$addToSet', [3, 4, 2, 5])])
        self.assertEqual(d['a'], [1, 2, 3, 4, 5])


class RenameTestCase(BaseTestCase):
    def test_renames(self):
        d = {'a': {'b': 'b'}}
        update_helper.update(d, [('a', '$rename', ('b', 'c'))])
        self.assertEqual(d['a'], {'c': 'b'})


class CallTestCase(BaseTestCase):
    def test_calls_for_dict(self):
        d = {'a': {'b': 'b'}}
        update_helper.update(d, [('a.setdefault', '$call', ('c', 'c'))])
        self.assertEqual(d['a'], {'b': 'b', 'c': 'c'})

    def test_calls_for_obj(self):
        class MyClass(object):
            def set_args(self, *args, **kwargs):
                self.args = args
                self.kwargs = kwargs

        o = MyClass()
        args = ['arg_1', 'arg_2']
        kwargs = {'k1': 'v1', 'k2': 'v2'}
        update_helper.update(o, [('set_args', '$call', args, kwargs)])
        self.assertEqual(list(o.args), list(args))
        self.assertEqual(o.kwargs, kwargs)


if __name__ == '__main__':
    unittest.main()
