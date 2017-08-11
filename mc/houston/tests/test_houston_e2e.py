import unittest


from mc.houston.houston import Houston
from mc.utils import test_utils as _mc_test_utils

from . import utils as _utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = Houston(cfg=self._generate_cfg())

    def _generate_cfg(self): return _utils.generate_test_cfg()

    def capture(self): return _mc_test_utils.capture()


class SanityCheckCommandTestCase(BaseTestCase):
    def test_calls_sanity_check_subcommand(self):
        with self.capture():
            self.houston.call_command('sanity_check')
