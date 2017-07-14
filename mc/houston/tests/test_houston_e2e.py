from contextlib import contextmanager
from io import StringIO
import sys
import unittest
from unittest.mock import MagicMock

from ..houston import Houston

@contextmanager
def capture():
  orig_stdout = sys.stdout
  try:
      sys.stdout = StringIO()
      sys.stdout.seek(0)
      yield sys.stdout.read()
  finally: sys.stdout = orig_stdout

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = Houston(cfg=MagicMock())

class SanityCheckCommandTestCase(BaseTestCase):
    def test_calls_sanity_check_subcommand(self):
        with capture(): self.houston.call_command('sanity_check')
