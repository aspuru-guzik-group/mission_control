import unittest

from mc.houston.tests import utils as _houston_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = _houston_test_utils.generate_test_houston()

    def _generate_cfg(self):
        return _houston_test_utils.generate_test_cfg()

    def _create_flow(self, flow_spec=None):
        flow_spec = flow_spec or {}
        return self.houston.run_command('create_flow', flow_spec=flow_spec)


class DefaultTestCase(BaseTestCase):
    def test_returns_summaries(self):
        self.assertEqual(
            self.houston.run_command('info'),
            {
                'flow': {'count': 0},
                'job': {'count': 0}
            }
        )
        self._create_flow()
        self.assertEqual(
            self.houston.run_command('info')['flow']['count'],
            1
        )


class KeyTestCase(BaseTestCase):
    def test_returns_record_for_key(self):
        flow = self._create_flow()
        self.assertEqual(
            self.houston.run_command('info', key=flow['key']),
            flow
        )
