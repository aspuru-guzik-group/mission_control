import collections
import unittest
from unittest.mock import call, MagicMock, patch

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
        self.submission_runner.create_workdir = MagicMock()
        self.submission_runner.run_workdir = MagicMock()
        self.submission_runner.move_to_outputs = MagicMock()
        self.expected_workdir_meta = \
                self.submission_runner.create_workdir.return_value
        self.submission_runner.run_submission()

    def test_creates_workdir(self):
        self.assertEqual(self.submission_runner.create_workdir.call_args,
                         call())

    def test_runs_workdir(self):
        self.assertEqual(self.submission_runner.run_workdir.call_args,
                         call(workdir_meta=self.expected_workdir_meta))

    def test_moves_workdir_to_outputs(self):
        self.assertEqual(
            self.submission_runner.move_to_outputs.call_args,
            call(src=self.expected_workdir_meta['dir'],
                 outputs_key=(
                     submission_runner.confgen_constants.CONFGEN_OUTPUTS_KEY)
                ))

class CreateWorkdirTestCase(BaseTestCase):
    @patch.object(submission_runner, 'WorkdirBuilder')
    def test_dispatches_to_build_workdir(self, MockWorkdirBuilder):
        self.submission_runner.generate_tmp_dir = MagicMock()
        workdir_meta = self.submission_runner.create_workdir()
        self.assertEqual(
            MockWorkdirBuilder.call_args,
            call(workdir=self.submission_runner.generate_tmp_dir.return_value,
                 workdir_params=self.submission_runner.get_workdir_params())
        )
        self.assertEqual(
            workdir_meta,
            MockWorkdirBuilder.return_value.build_workdir.return_value)

class RunWorkdirTestCase(BaseTestCase):
    @patch.object(submission_runner, 'os')
    @patch.object(submission_runner, 'subprocess')
    def test_runs_workdir_entrypoint(self, mock_subprocess, mock_os):
        workdir_meta = MagicMock()
        self.submission_runner.run_workdir(workdir_meta=workdir_meta)
        expected_entrypoint_path = mock_os.path.join(workdir_meta['dir'],
                                                     workdir_meta['entrypoint'])
        expected_cmd = ['/bin/bash', expected_entrypoint_path]
        expected_env = {
            **mock_os.environ,
            **(self.submission_runner.cfg.get(
                'a2g2.jobs.confgen.confgen', {}).get('env_vars', {}))
        }
        self.assertEqual(mock_subprocess.run.call_args,
                         call(expected_cmd, check=True, env=expected_env))
