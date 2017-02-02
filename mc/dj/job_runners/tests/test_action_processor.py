import unittest
from unittest.mock import call, MagicMock, Mock

from ..action_processor import ActionProcessor, ActionCtx


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.action_processor = ActionProcessor()

class InitTestCase(BaseTestCase):
    def setUp(self):
        self.handlers = {i: Mock() for i in range(3)}
        self.action_processor = ActionProcessor(handlers=self.handlers)

    def test_registers_handlers(self):
        registered_handlers = set(self.action_processor.handlers.items())
        self.assertTrue(registered_handlers.issuperset(self.handlers.items()))

    def test_has_default_set_ctx_value_handler(self):
        self.assertEqual(self.action_processor.handlers['set_ctx_value'],
                         self.action_processor.set_ctx_value_handler)

class SetCtxValueHandler(BaseTestCase):
    def test_wraps_set_ctx_value(self):
        ctx = MagicMock()
        self.action_processor.set_ctx_value = Mock()
        params = {str(i): Mock() for i in range(3)}
        self.action_processor.set_ctx_value_handler(params=params, ctx=ctx)
        self.assertEqual(self.action_processor.set_ctx_value.call_args,
                         call(ctx=ctx, **params))

class SetCtxValueTestCase(BaseTestCase):
    def test_calls_ctx_set_target_with_transformed_value(self):
        ctx = Mock()
        target = Mock()
        value = Mock()
        self.action_processor.set_ctx_value(ctx=ctx, target=target, value=value)
        self.assertEqual(ctx.transform_value.call_args, call(value=value))
        self.assertEqual(
            ctx.set.call_args,
            call(target=target, value=ctx.transform_value.return_value)
        )

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

class ActionCtxBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.ctx_to_wrap = Mock()
        self.ctx = ActionCtx(ctx_to_wrap=self.ctx_to_wrap)

class TransformValue(ActionCtxBaseTestCase):
    def setUp(self):
        super().setUp()
        self.ctx.render_template = Mock()

    def test_renders_template_values(self):
        value = {'template': 'some template'}
        result = self.ctx.transform_value(value)
        self.assertEqual(self.ctx.render_template.call_args,
                         call(value['template']))
        self.assertEqual(result, self.ctx.render_template.return_value)

    def test_does_not_render_non_template_values(self):
        value = 'some non-template value'
        result = self.ctx.transform_value(value)
        self.assertEqual(self.ctx.render_template.call_count, 0)
        self.assertEqual(result, value)
