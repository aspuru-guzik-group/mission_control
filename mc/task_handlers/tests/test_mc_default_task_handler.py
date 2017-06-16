import collections
import unittest
from unittest.mock import call, patch, MagicMock

from .. import mc_default_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = collections.defaultdict(MagicMock)
        self.task_ctx = collections.defaultdict(MagicMock, {'task': self.task})
        self.task_handler = mc_default_task_handler.McDefaultTaskHandler(
            task_ctx=self.task_ctx)
        self.args = [MagicMock() for i in range(3)]
        self.kwargs = {'kwargs_%s' % i: MagicMock() for i in range(3)}

    def setup_mocks(self, attrs=None):
        patchers = {attr: patch.object(mc_default_task_handler, attr)
                    for attr in attrs}
        mocks = {}
        for attr, patcher in patchers.items():
            mocks[attr] = patcher.start()
            self.addCleanup(patcher.stop)
        return mocks

    def setup_handler_mocks(self, attrs=None):
        for attr in attrs: setattr(self.task_handler, attr, MagicMock())


class DefaultTaskTypeToHandlerDotSpecTestCase(BaseTestCase):
    def _assert_task_type_to_dot_spec(self, task_type=None, expected_dot_spec=None):
        result = self.task_handler.default_task_type_to_handler_dot_spec(
            task_type=task_type)
        self.assertEqual(result, expected_dot_spec)

    def test_ignores_specs_w_colon(self):
        task_type = 'some.spec.with:a_colon'
        self._assert_task_type_to_dot_spec(task_type=task_type,
                                           expected_dot_spec=task_type)

    def test_replaces_task_w_task_handler(self):
        task_type = 'some.task'
        self._assert_task_type_to_dot_spec(
            task_type=task_type,
            expected_dot_spec='some.task_handler:TaskHandler'
        )

    def test_appends_task_handler(self):
        task_type = 'some.task.foo'
        self._assert_task_type_to_dot_spec(
            task_type=task_type,
            expected_dot_spec='some.task_handler.foo_task_handler:TaskHandler'
        )

    def test_appends_colon_TaskHandler(self):
        task_type = 'some.task.foo_task_handler'
        self._assert_task_type_to_dot_spec(
            task_type=task_type,
            expected_dot_spec='some.task_handler.foo_task_handler:TaskHandler'
        )

class _TickTaskTestCase(BaseTestCase): pass

class InterpolateAndTickTaskCtxTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.task = MagicMock()
        self.setup_handler_mocks(attrs=['get_interpolated_task_ctx',
                                        'tick_task_ctx'])
        self.task_handler.interpolate_and_tick_task_ctx(
            *self.args, **self.kwargs)
        self.expected_interpolated_task_ctx = (
            self.task_handler.get_interpolated_task_ctx.return_value)

    def test_gets_interpolated_task_ctx(self):
        self.assertEqual(self.task_handler.get_interpolated_task_ctx.call_args,
                         call())

    def test_calls_tick_task_ctx(self):
        self.assertEqual(
            self.task_handler.tick_task_ctx.call_args,
            call(*self.args, task_ctx=self.expected_interpolated_task_ctx,
                 **self.kwargs)
        )

    def test_propagates_task_ctx_changes(self):
        self.assertEqual(self.task_handler.task.update.call_args,
                         call(self.expected_interpolated_task_ctx['task']))

class TickTaskCtxTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_handler_mocks(attrs=['get_handler_for_task_ctx'])
        self.task_handler.tick_task_ctx(
            *self.args, task_ctx=self.task_ctx, **self.kwargs)

    def test_gets_handler_for_task_ctx(self):
        self.assertEqual(
            self.task_handler.get_handler_for_task_ctx.call_args,
            call(task_ctx=(self.task_ctx))
        )

    def test_calls_handler_tick_task(self):
        expected_handler = (self.task_handler.get_handler_for_task_ctx
                            .return_value)
        self.assertEqual(
            expected_handler.tick_task.call_args,
            call(*self.args, task_ctx=self.task_ctx, **self.kwargs)
        )

class GetHandlerForTaskCtxTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task = collections.defaultdict(MagicMock)

    def _get_handler(self, task=None):
        return self.task_handler.get_handler_for_task_ctx(
            task_ctx={**self.task_ctx, 'task': task})

    def assert_task_type_to_handler(self, task_type=None,
                                    expected_handler=None):
        handler = self._get_handler(task={'task_type': task_type})
        self.assertEqual(handler, expected_handler)

    def test_log(self):
        self.assert_task_type_to_handler(
            'log', mc_default_task_handler.LogTaskHandler)

    def test_noop(self):
        self.assert_task_type_to_handler(
            'noop', mc_default_task_handler.NoOpTaskHandler)

    def test_print(self):
        self.assert_task_type_to_handler(
            'print', mc_default_task_handler.PrintTaskHandler)

    def test_wire(self):
        self.assert_task_type_to_handler(
            'wire', mc_default_task_handler.WireTaskHandler)

    def test_dot_spec(self):
        self.setup_handler_mocks(attrs=['task_type_to_handler_dot_spec',
                                        'load_from_dot_spec'])
        handler = self._get_handler(task=self.task)
        self.assertEqual(
            self.task_handler.task_type_to_handler_dot_spec.call_args,
            call(self.task['task_type'])
        )
        self.assertEqual(
            self.task_handler.load_from_dot_spec.call_args,
            call(dot_spec=(self.task_handler.task_type_to_handler_dot_spec
                           .return_value))
        )
        self.assertEqual(handler,
                         self.task_handler.load_from_dot_spec.return_value)

class GetInterpolatedTaskCtxTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_ctx.update({'key_%s' % i: MagicMock() for i in range(3)})
        self.setup_handler_mocks(attrs=['get_interpolated_task'])
        self.result = self.task_handler.get_interpolated_task_ctx()

    def test_returns_task_ctx_w_interpolated_task(self):
        expected = {
            **self.task_ctx,
            'task': self.task_handler.get_interpolated_task.return_value
        }
        self.assertEqual(self.task_handler.get_interpolated_task.call_args,
                         call())
        self.assertEqual(self.result, expected)

class GetInterpolatedTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mocks = self.setup_mocks(attrs=['copy'])
        self.setup_handler_mocks(attrs=[
            'traverse_obj', 'value_has_interpolations', 'interpolate_value'])
        self.mock_parent = {}
        self.mock_traversal_items = []
        for i in range(4):
            key = 'key_%s' % i
            value = {'id': i, 'interpolate': (i % 2) == 0}
            self.mock_traversal_items.append((self.mock_parent, key, value))
        self.task_handler.traverse_obj.return_value = self.mock_traversal_items
        def mock_has_interpolations(value=None): return value['interpolate']
        self.task_handler.value_has_interpolations = mock_has_interpolations
        self.result = self.task_handler.get_interpolated_task()

    def test_copies_task(self):
        self.assertEqual(self.mocks['copy'].deepcopy.call_args,
                         call(self.task))

    def test_traverses_copy(self):
        self.assertEqual(self.task_handler.traverse_obj.call_args,
                         call(self.mocks['copy'].deepcopy.return_value))

    def test_interpolates_values_that_have_interpolations(self):
        expected_interpolations = {
            key: self.task_handler.interpolate_value.return_value
            for parent, key, value in self.mock_traversal_items
            if value['interpolate']
        }
        self.assertEqual(self.mock_parent, expected_interpolations)

    def test_returns_copy(self):
        self.assertEqual(self.result, self.mocks['copy'].deepcopy.return_value)

class ValueHasInterpolationsTestCase(BaseTestCase):
    def _assert_has_interpolations(self, value=None, expected=None):
        self.assertEqual(self.task_handler.value_has_interpolations(value),
                         expected)

    def test_true_for_money_string(self):
        self._assert_has_interpolations('$money!', True)

    def test_true_for_escaped_money_string(self):
        self._assert_has_interpolations('\$money!', True)

    def test_false_for_non_money(self):
        self._assert_has_interpolations('\no dinero baby :(', False)

    def test_false_for_non_string(self):
        self._assert_has_interpolations(MagicMock(), False)

class InterpolateValueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.ctx = MagicMock()
        self.setup_handler_mocks(attrs=['get_ctx_value'])

    def _interpolate(self, value=None):
        return self.task_handler.interpolate_value(value=value, ctx=self.ctx)

    def _assert_interpolated(self, value=None, expected_interpolated=None):
        self.assertEqual(self._interpolate(value), expected_interpolated)

    def test_money_value(self):
        value = '$money'
        self._assert_interpolated(
            value=value,
            expected_interpolated=self.task_handler.get_ctx_value.return_value
        )
        self.assertEqual(self.task_handler.get_ctx_value.call_args,
                         call(ctx=self.ctx, dot_spec=value[1:]))

    def test_escaped_money_value(self):
        value = '\$escaped_money'
        self._assert_interpolated(value=value, expected_interpolated=value[1:])

if __name__ == '__main__':
    unittest.main()
