from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import job_engine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job_engine = job_engine.JobEngine()
        self.job = MagicMock()
        self.cfg = defaultdict(MagicMock)
        self.output_dir = 'some_output_dir'

class ExecuteCommandTestCase(BaseTestCase):
    def test_build_submission_command(self):
        self.assert_handles_command(command='build_submission')

    def assert_handles_command(self, command=None):
        mock_handler = MagicMock()
        setattr(self.job_engine, command, mock_handler)
        kwargs = {'job': self.job, 'cfg': self.cfg, 'extra': 'extra'}
        result = self.job_engine.execute_command(command=command, **kwargs)
        self.assertEqual(mock_handler.call_args, call(**kwargs))
        self.assertEqual(result, mock_handler.return_value)

    def test_run_submission_command(self):
        self.assert_handles_command(command='run_submission')

@patch.object(job_engine, 'os')
class BuildSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.output_dir = MagicMock()
        self.job_engine.get_job_module = MagicMock()
        self.job_engine.write_submission_meta = MagicMock()
        self.expected_module = self.job_engine.get_job_module.return_value

    def _build_submission(self):
        return self.job_engine.build_submission(
            job=self.job, cfg=self.cfg, output_dir=self.output_dir)

    def test_dispatches_to_job_module(self, *args):
        self._build_submission()
        self.assertEqual(
            self.job_engine.get_job_module.call_args,
            call(job=self.job, cfg=self.cfg)
        )
        self.assertEqual(
            self.expected_module.build_submission.call_args,
            call(job=self.job, cfg=self.cfg, output_dir=self.output_dir)
        )

    def test_writes_and_returns_submission_meta(self, mock_os):
        result = self._build_submission()
        expected_submission_meta = {
            'job': self.job, 'cfg': self.cfg, 'dir': self.output_dir,
            **self.expected_module.build_submission.return_value
        }
        self.assertEqual(
            self.job_engine.write_submission_meta.call_args,
            call(submission_meta=expected_submission_meta, _dir=self.output_dir)
        )
        self.assertEqual(result, expected_submission_meta)

class GetJobModuleCase(BaseTestCase):
    @patch.object(job_engine, 'importlib')
    def test_imports_expected_module(self, patched_importlib):
        self.job_engine.get_job_module_name = MagicMock()
        result = self.job_engine.get_job_module(job=self.job, cfg=self.cfg)
        self.assertEqual(
            patched_importlib.import_module.call_args,
            call(self.job_engine.get_job_module_name.return_value)
        )
        self.assertEqual(result, patched_importlib.import_module.return_value)

class RunSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_dir = MagicMock()
        self.submission_meta = MagicMock()
        self.job_engine.read_submission_meta = \
                MagicMock(return_value=self.submission_meta)
        self.job_module = MagicMock()
        self.job_engine.get_job_module = MagicMock(return_value=self.job_module)
        self.job_engine.run_submission(submission_dir=self.submission_dir)

    def test_reads_submission_meta(self):
        self.assertEqual(self.job_engine.read_submission_meta.call_args,
                         call(_dir=self.submission_dir))

    def test_dispatches_to_job_module(self):
        self.assertEqual(
            self.job_engine.get_job_module.call_args,
            call(job=self.submission_meta['job'],
                 cfg=self.submission_meta['cfg'])
        )
        self.assertEqual(self.job_module.run_submission.call_args,
                         call(submission_meta=self.submission_meta))
