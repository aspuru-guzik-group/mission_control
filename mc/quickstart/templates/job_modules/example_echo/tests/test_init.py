import unittest

from ... import confgen

class BuildSubmissionTestCase(unittest.TestCase):
    def test_has_build_job_submission(self, *args):
        self.assertTrue(confgen.build_job_submission is not None)

class RunSubmissionTestCase(unittest.TestCase):
    def test_has_run_job_submission(self, *args):
        self.assertTrue(confgen.run_job_submission is not None)
