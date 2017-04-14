import unittest

from ... import query


class BuildSubmissionTestCase(unittest.TestCase):
    def test_aliases_odyssey_submission_builder(self, *args):
        self.assertEqual(
            query.build_job_submission,
            query.odyssey_submission_builder.build_job_submission
        )

class RunSubmissionTestCase(unittest.TestCase):
    def test_aliases_submission_runner(self, *args):
        self.assertEqual(query.run_job_submission,
                         query.query_submission_runner.run_job_submission)
