import unittest

from ... import flow 


class BuildSubmissionTestCase(unittest.TestCase):
    def test_aliases_odyssey_submission_builder(self, *args):
        self.assertEqual(
            flow.build_job_submission,
            flow.odyssey_submission_builder.build_job_submission
        )

class RunSubmissionTestCase(unittest.TestCase):
    def test_aliases_submission_runner(self, *args):
        self.assertEqual(flow.run_job_submission,
                         flow.flow_submission_runner.run_job_submission)
