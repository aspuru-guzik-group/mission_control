from collections import defaultdict
import json
import os
import textwrap
import tempfile
import unittest
from unittest.mock import call, patch, MagicMock

from .. import submission_builder


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = defaultdict(MagicMock)
        self.cfg = defaultdict(MagicMock)
        self.submission_dir = MagicMock()
        self.submission_builder = submission_builder.SubmissionBuilder(
            job=self.job, cfg=self.cfg, submission_dir=self.submission_dir)

@patch.object(submission_builder, 'OdysseyJobDirBuilder')
class BuildJobSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_builder.ensure_dir = MagicMock()
        self.submission_builder.generate_entrypoint_body = MagicMock()

    def test_dispatches_to_odyssey_builder(self, MockOdysseyBuilder):
        self.submission_builder.build_submission()
        self.assertEqual(
            MockOdysseyBuilder.build_dir.call_args,
            call(dir_spec={'entrypoint_body': \
                           self.submission_builder.generate_entrypoint_body()},
                 output_dir=self.submission_dir,
                 submission_meta_file_name=(
                     self.submission_builder.submission_meta_file_name)
                )
        )

    def test_returns_submission_meta(self, MockOdysseyBuilder):
        submission_meta = self.submission_builder.build_submission()
        expected_submission_meta = MockOdysseyBuilder.build_dir.return_value
        self.assertEqual(submission_meta, expected_submission_meta)

class GenerateEntrypointBodyTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_builder.write_json_params = MagicMock(return_value={})

    def test_generates_expected_entrypoint_body(self):
        job_engine_cfg = self.submission_builder.cfg['job_engine']
        expected_entrypoint_body = textwrap.dedent(
            """
            {job_engine_preamble}
            python -m {job_engine_module} {job_engine_command} {cli_params}
            """
        ).strip().format(
            job_engine_preamble=job_engine_cfg.get('entrypoint_preamble', ''),
            job_engine_module=job_engine_cfg['engine_module'],
            job_engine_command='run_job_submission',
            cli_params=self.submission_builder.params_to_cli_args(
                params={
                    **(self.submission_builder.write_json_params()),
                    'submission': os.path.join(
                        self.submission_builder.submission_dir,
                        self.submission_builder.submission_meta_file_name
                    )
                }
            )
        )
        self.assertEqual(self.submission_builder.generate_entrypoint_body(),
                         expected_entrypoint_body)

class WriteJsonParamsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_builder.cfg = {'cfg_key': 'cfg_value'}
        self.submission_builder.job = {'job_key': 'job_value'}
        self.submission_builder.submission_dir = tempfile.mkdtemp()

    def test_writes_expected_json_files(self):
        actual_params = self.submission_builder.write_json_params()
        actual_contents = {}
        expected_contents = {}
        expected_params = {}
        for param in ['job', 'cfg']:
            json_filename = param + '.json'
            expected_params[param] = json_filename
            json_path = os.path.join(
                self.submission_builder.submission_dir, json_filename)
            actual_contents[param] = open(json_path).read()
            expected_contents[param] = json.dumps(
                getattr(self.submission_builder, param))
        self.assertEqual(actual_contents, expected_contents)
        self.assertEqual(actual_params, expected_params)
