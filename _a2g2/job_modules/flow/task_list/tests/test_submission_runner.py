import collections
import unittest
from unittest.mock import MagicMock

from .. import submission_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = MagicMock()
        self.cfg = MagicMock()
        self.submission = collections.defaultdict(MagicMock)
        self.submission_runner = submission_runner.SubmissionRunner(
            job=self.job, cfg=self.cfg, submission=self.submission)

class GetFlowSpecTestCase(BaseTestCase):
    def test_returns_expected_flow_spec(self):
        result = self.submission_runner.get_flow_spec()
        job_params = self.job['job_spec']['job_params']
        expected_flow_spec = {
            'data': job_params.get('flow_data'),
            'node_specs': [
                {'node': {'node_tasks': job_params['tasks']},
                 'precursor_keys': ['ROOT']}
            ]
        }
        self.assertEqual(result, expected_flow_spec)
