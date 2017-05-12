import unittest

from ... import confgen


class BuildSubmissionTestCase(unittest.TestCase):
    def test_aliases_odyssey_submission_builder(self, *args):
        self.assertEqual(
            confgen.build_job_submission,
            confgen.odyssey_submission_builder.build_job_submission
        )

class RunSubmissionTestCase(unittest.TestCase):
    def test_aliases_submission_runner(self, *args):
        self.assertEqual(confgen.run_job_submission,
                         confgen.confgen_submission_runner.run_job_submission)
