import unittest
from unittest.mock import MagicMock

from .. import houston


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.subcommands_registry = {
            'subcommand_%s' % i: MagicMock(return_value=None)
            for i in range(3)
        }
        self.houston = houston.HoustonCommand(
            subcommands_registry=self.subcommands_registry)

class RunTestCase(BaseTestCase):
    def test_dispatches_to_subcommand(self):
        subcommand = list(self.subcommands_registry.keys())[0]
        self.houston.run_from_argv(argv=['houston', subcommand])
        expected_subcommand_fn = self.subcommands_registry[subcommand]
        call_kwargs = expected_subcommand_fn.call_args[1]
        self.assertEqual(call_kwargs['get_cfg'], self.houston._get_cfg)
