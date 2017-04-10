import unittest

from ... import qchem 


class BuildSubmissionTestCase(unittest.TestCase):
    def test_aliases_odyssey_submission_builder(self, *args):
        self.assertEqual(
            qchem.build_job_submission,
            qchem.odyssey_submission_builder.build_job_submission
        )

class RunSubmissionTestCase(unittest.TestCase):
    def test_aliases_submission_runner(self, *args):
        self.assertEqual(qchem.run_job_submission,
                         qchem.qchem_submission_runner.run_job_submission)
