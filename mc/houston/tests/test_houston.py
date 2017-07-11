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

    def mockify_houston_attrs(self, attrs=None):
        for attr in attrs: setattr(self.houston, attr, MagicMock())

class RunTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_houston_attrs(attrs=['_get_load_cfg_fn'])
        self.subcommand = list(self.subcommands_registry.keys())[0]
        self.houston.run_from_argv(argv=['houston', self.subcommand])

    def test_dispatches_to_subcommand(self):
        expected_subcommand_fn = self.subcommands_registry[self.subcommand]
        call_kwargs = expected_subcommand_fn.call_args[1]
        self.assertEqual(call_kwargs['load_cfg'],
                         self.houston._get_load_cfg_fn.return_value)
