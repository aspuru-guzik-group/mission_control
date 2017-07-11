import unittest
from unittest import call, MagicMock

from .. import _utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.cfg = MagicMock()
        self.utils = _utils.HoustonSubcomandUtils(
            get_cfg=MagicMock(return_value=self.cfg)
        )

    def mockify_utils_attrs(self, attrs=None):
        for attr in attrs: setattr(self.utils, attr, MagicMock())

class _GetFlowRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_get_mc_client',
                                             '_get_mc_queues'])
        self.mod_mocks = self.mockify_module_attrs(attrs=['FlowRunner'])
        self.result = self.subcommand._get_flow_runner()

    def test_gets_mc_client(self):
        self.assertEqual(self.subcommand._get_mc_client.call_args, call())
        
    def test_gets_flow_queue_key(self):
        self.assertEqual(self.subcommand._get_mc_queues.call_args, call())

    def test_constructs_and_returns_flow_runner(self):
        self.assertEqual(self.mod_mocks['FlowRunner'].call_args,
                         call())
        self.assertEqual(self.result, self.mod_mocks['FlowRunner'].return_value)

class _GetMcClientTestCase(BaseTestCase):
    def test_something(self): self.fail()

class _GetMcQueuesTestCase(BaseTestCase):
    def test_something(self): self.fail()
