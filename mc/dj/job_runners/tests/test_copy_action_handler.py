import unittest
from unittest.mock import call, Mock

from ..copy_action_handler import CopyActionHandler


class CopyActionHandlerTestCase(unittest.TestCase):
    def setUp(self):
        self.transfer_client = Mock()
        self.mkdtemp = Mock(return_value='/mock_tmp_dir')
        self.handler = CopyActionHandler(transfer_client=self.transfer_client,
                                         mkdtemp=self.mkdtemp)
        self.ctx = Mock()
        self.params = {'src': Mock(), 'dest': Mock()}
        self.output = self.handler(params=self.params, ctx=self.ctx)
        self.expected_intermediate_dir = self.handler.get_intermediate_dir()

    def test_copies_src_to_tmp_dir(self):
        self.assertEqual(
            self.transfer_client.copy.call_args_list[0],
            call(src=self.params['src'], dest=self.expected_intermediate_dir)
        )

    def test_copies_tmp_dir_to_dest(self):
        self.assertEqual(
            self.transfer_client.copy.call_args_list[1],
            call(src=self.expected_intermediate_dir, dest=self.params['dest'])
        )

    def test_returns_transfer_output(self):
        self.assertEqual(self.output, self.transfer_client.copy.return_value)
