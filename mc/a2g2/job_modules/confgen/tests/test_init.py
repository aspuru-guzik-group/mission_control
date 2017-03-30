import unittest
from unittest.mock import call, MagicMock, patch

from ... import confgen


class BuildSubmissionTestCase(unittest.TestCase):
    @patch.object(confgen, 'get_submission_builder')
    def test_dispatches_to_submission_builder(self, *args):
        job = MagicMock()
        cfg = MagicMock()
        submission_dir = MagicMock()
        submission_meta = confgen.build_job_submission(
            job=job, cfg=cfg, submission_dir=submission_dir)
        self.assertEqual(
            confgen.get_submission_builder.call_args,
            call(job=job, cfg=cfg, submission_dir=submission_dir)
        )
        expected_builder = confgen.get_submission_builder.return_value
        expected_build_submission_fn = expected_builder.build_submission
        self.assertEqual(expected_build_submission_fn.call_args, call())
        self.assertEqual(submission_meta,
                         expected_build_submission_fn.return_value)

class RunSubmissionTestCase(unittest.TestCase):
    @patch.object(confgen, 'get_submission_runner')
    def test_dispatches_to_submission_runner(self, *args):
        submission = MagicMock()
        confgen.run_job_submission(submission=submission)
        self.assertEqual(confgen.get_submission_runner.call_args,
                         call(submission=submission))
        expected_runner = confgen.get_submission_runner.return_value
        self.assertEqual(expected_runner.run_submission.call_args, call())
