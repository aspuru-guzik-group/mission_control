import json
import unittest

from mc.houston.houston import Houston
from mc.houston import test_utils as _houston_test_utils
from mc.utils import test_utils as _mc_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = Houston(cfg=self._generate_cfg())
        self.houston.utils.ensure_db()

    def capture(self): return _mc_test_utils.capture()

    def _generate_cfg(self):
        return _houston_test_utils.generate_test_cfg()

    def _create_flow(self, flow_spec=None):
        flow_spec = flow_spec or {}
        with self.capture() as stdout:
            self.houston.call_command(
                'create_flow', flow_spec=json.dumps(flow_spec))
        return json.loads(stdout.read())

    def _capture_info(self, *args, **kwargs):
        with self.capture() as stdout:
            self.houston.call_command('info', *args, **kwargs)
        parsed_info = json.loads(stdout.read())
        return parsed_info

class DefaultTestCase(BaseTestCase):
    def test_returns_summaries(self):
        self.assertEqual(self._capture_info(),
            {
                'flow': {'count': 0},
                'job': {'count': 0}
            }
        )
        self._create_flow()
        self.assertEqual(self._capture_info()['flow']['count'], 1)

class KeyTestCase(BaseTestCase):
    def test_returns_record_for_key(self):
        flow = self._create_flow()
        self.assertEqual(self._capture_info(key=flow['key']), flow)
