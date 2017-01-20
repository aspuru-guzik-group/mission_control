import unittest
from unittest.mock import call, MagicMock

from ..copy_action_handler import CopyActionHandler

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.transfer_client = MagicMock()
        self.render_tpl = MagicMock()
        self.handler = CopyActionHandler(transfer_client=self.transfer_client,
                                         render_tpl=self.render_tpl)
        self.ctx = MagicMock()
        self.params = {'src': MagicMock(), 'dest': MagicMock()}

class CopyActionHandlerTestCase(BaseTestCase):
    def test_calls_transfer_client_with_formatted_params(self):
        self.handler.format_params = MagicMock()
        self.handler(params=self.params, ctx=self.ctx)
        expected_formatted_params = self.handler.format_params.return_value
        expected_copy_params = {key: expected_formatted_params[key]
                                for key in ['src', 'dest']}
        self.assertEqual(self.transfer_client.copy.call_args,
                         call(**expected_copy_params))

class FormatParamsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.formatted_params = self.handler.format_params(params=self.params,
                                                           ctx=self.ctx)
    def test_renders_param_templates_with_ctx(self):
        expected_call_args_list = [call(tpl=self.params[key], ctx=self.ctx)
                                   for key in ['src', 'dest']]
        self.assertEqual(self.render_tpl.call_args_list,
                         expected_call_args_list)

    def test_returns_rendered_templates(self):
        expected_formatted_params = {key: self.render_tpl.return_value
                                     for key in ['src', 'dest']}
        self.assertEqual(self.formatted_params, expected_formatted_params)

