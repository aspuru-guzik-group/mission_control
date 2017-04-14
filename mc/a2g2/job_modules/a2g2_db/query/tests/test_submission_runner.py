import collections
import unittest
from unittest.mock import call, patch, MagicMock

from .. import submission_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = MagicMock()
        self.cfg = MagicMock()
        self.submission = collections.defaultdict(MagicMock)
        self.submission_runner = submission_runner.SubmissionRunner(
            job=self.job, cfg=self.cfg, submission=self.submission)

class RunSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_runner.execute_query = MagicMock()
        self.submission_runner.write_query_results = MagicMock()
        self.submission_runner.run_submission()

    def test_executes_query(self):
        expected_query_params = \
                self.job['job_spec']['job_params']['query_params']
        self.assertEqual(self.submission_runner.execute_query.call_args,
                         call(query_params=expected_query_params))

    def test_writes_query_results(self):
        expected_query_results = \
                self.submission_runner.execute_query.return_value
        self.assertEqual(self.submission_runner.write_query_results.call_args,
                         call(query_results=expected_query_results))

@patch.object(submission_runner.a2g2_db_utils, 'generate_a2g2_client')
class ExecuteQueryTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query_params = MagicMock()

    def test_dispatches_to_a2g2_client(self, mock_generate_a2g2_client):
        query_results = self.submission_runner.execute_query(
            query_params=self.query_params)
        self.assertEqual(
            submission_runner.a2g2_db_utils.generate_a2g2_client.call_args,
            call(cfg=self.submission_runner.cfg)
        )
        expected_client = submission_runner.a2g2_db_utils.generate_a2g2_client\
                .return_value
        self.assertEqual(expected_client.execute_query.call_args,
                         call(**self.query_params))
        self.assertEqual(query_results,
                         expected_client.execute_query.return_value)

@patch.object(submission_runner, 'sys')
@patch.object(submission_runner, 'json')
class WriteQueryResultsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query_results = MagicMock()

    def test_writes_jsonified_query_results_to_stdout(self, mock_json,
                                                      mock_sys):
        self.submission_runner.write_query_results(
            query_results=self.query_results)
        self.assertEqual(mock_json.dumps.call_args, call(self.query_results))
        expected_serialization = mock_json.dumps.return_value
        self.assertEqual(mock_sys.stdout.write.call_args,
                         call(expected_serialization))
