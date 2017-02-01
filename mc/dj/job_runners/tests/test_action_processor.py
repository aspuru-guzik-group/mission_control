import unittest
from unittest.mock import call, MagicMock, Mock

from ..action_processor import ActionProcessor


class BaseTestCase(unittest.TestCase):
    pass

class InitTestCase(BaseTestCase):
    def test_registers_handlers(self):
        handlers = {i: Mock() for i in range(3)}
        action_processor = ActionProcessor(handlers=handlers)
        contains_handlers = set(action_processor.handlers.items()).issuperset(
            set(handlers.items()))
        self.assertTrue(contains_handlers)

    def test_has_default_set_ctx_value_handler(self):
        action_processor = ActionProcessor()
        self.assertEqual(action_processor.handlers['set_ctx_value'],
                         action_processor.set_ctx_value_handler)

class SetCtxValueHandler(BaseTestCase):
    def test_wraps_set_ctx_value(self):
        ctx = MagicMock()
        action_processor = ActionProcessor()
        action_processor.set_ctx_value = Mock()
        params = {str(i): Mock() for i in range(3)}
        action_processor.set_ctx_value_handler(params=params, ctx=ctx)
        self.assertEqual(action_processor.set_ctx_value.call_args,
                         call(ctx=ctx, **params))

class SetCtxValueTestCase(BaseTestCase):
    def test_calls_ctx_set(self):
        ctx = MagicMock()
        target = Mock()
        value = Mock()
        action_processor = ActionProcessor()
        action_processor.set_ctx_value(ctx=ctx, target=target, value=value)
        self.assertEqual(ctx.set.call_args,
                         call(target=target, value=value))

class ProcessActionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.ctx = Mock()
        self.handler_key = 'handler_key'
        self.handler = Mock()
        self.handlers = {self.handler_key: self.handler}
        self.action_processor = ActionProcessor(handlers=self.handlers)
        self.action = {
            'action': self.handler_key,
            'params': Mock(),
            'output_to_ctx_target': 'some.ctx.target',
        }
        self.action_processor._wrap_ctx_for_action = MagicMock()

    def test_dispatches_to_handler_with_wrapped_ctx(self):
        self.action_processor.process_action(action=self.action, ctx=self.ctx)
        self.assertEqual(
            self.handler.call_args,
            call(params=self.action['params'],
                 ctx=self.action_processor._wrap_ctx_for_action.return_value)
        )

    def test_handles_store_output_to_ctx(self):
        self.action_processor.set_ctx_value = Mock()
        self.action_processor.process_action(action=self.action, ctx=self.ctx)
        self.assertEqual(
            self.action_processor.set_ctx_value.call_args,
            call(ctx=self.action_processor._wrap_ctx_for_action.return_value,
                 target=self.action['output_to_ctx_target'],
                 value=self.handler.return_value)
        )
