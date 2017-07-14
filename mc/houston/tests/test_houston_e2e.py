from contextlib import contextmanager
from io import StringIO
import json
import sys
import unittest

from jobman.engines.local_engine import LocalEngine

from ..houston import Houston


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = Houston(cfg=self._generate_cfg())
        self.houston.utils.ensure_db()

    @contextmanager
    def capture(self):
        orig_stdout = sys.stdout
        try:
            tmp_stdout = StringIO()
            sys.stdout = tmp_stdout
            yield tmp_stdout
            tmp_stdout.seek(0)
        finally: sys.stdout = orig_stdout

    def _generate_cfg(self):
        cfg = {
            'MC_DB_URI': 'sqlite://',
            'JOBMAN_CFG': {
                'jobman_db_uri': 'sqlite://',
                'engine': LocalEngine(db_uri='sqlite://')
            }
        }
        return cfg

class SanityCheckCommandTestCase(BaseTestCase):
    def test_calls_sanity_check_subcommand(self):
        with self.capture(): self.houston.call_command('sanity_check')

class InfoTestCase(BaseTestCase):
    def _create_flow(self, flow_spec=None):
        flow_spec = flow_spec or {}
        with self.capture():
            self.houston.call_command('create_flow',
                                      flow_spec=json.dumps(flow_spec))

    def _get_info(self):
        with self.capture() as stdout: self.houston.call_command('info')
        parsed_info = json.loads(stdout.read())
        return parsed_info

    def test_returns_general_info(self):
        self.assertEqual(self._get_info(),
            {
                'Flow': {'count': 0},
                'Job': {'count': 0}
            }
        )
        self._create_flow()
        self.assertEqual(self._get_info()['Flow']['count'], 1)
