import unittest
from unittest.mock import MagicMock

from .. import houston_command


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.subcommands = {
            'subcommand_%s' % i: MagicMock(return_value=None)
            for i in range(3)
        }
        self.houston = houston_command.HoustonCommand(
            subcommands=self.subcommands)

    def mockify_houston_attrs(self, attrs=None):
        for attr in attrs:
            setattr(self.houston, attr, MagicMock())


class RunTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_houston_attrs(attrs=['_load_cfg_from_path'])
        self.subcommand = list(self.subcommands.keys())[0]
        self.houston.run_from_argv(argv=['houston', self.subcommand])

    def test_dispatches_to_subcommand(self):
        expected_subcommand = self.subcommands[self.subcommand]
        self.assertEqual(len(expected_subcommand.call.call_args_list), 1)
