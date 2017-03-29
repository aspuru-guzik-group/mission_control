from collections import defaultdict
import json
import os
import textwrap
import tempfile
import time
import unittest
from unittest.mock import call, MagicMock

from .. import submission_builder


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = defaultdict(MagicMock)
        self.cfg = defaultdict(MagicMock)
        self.submission_dir = MagicMock()
        self.submission_builder = submission_builder.SubmissionBuilder(
            job=self.job, cfg=self.cfg, submission_dir=self.submission_dir)

class BuildSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_builder.ensure_dir = MagicMock()
        self.submission_builder.write_entrypoint = MagicMock()
        self.submission_builder.entrypoint_file_name = MagicMock()

    def test_writes_entrypoint(self):
        self.submission_builder.build_submission()
        self.assertEqual(self.submission_builder.write_entrypoint.call_args,
                         call())

    def test_returns_submission_meta(self):
        submission_meta = self.submission_builder.build_submission()
        expected_submission_meta = {
            'dir': self.submission_dir,
            'entrypoint': (
                self.submission_builder.write_entrypoint.return_value)
        }
        self.assertEqual(submission_meta, expected_submission_meta)

class WriteEntrypointTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.tmpdir = tempfile.mkdtemp()
        self.entrypoint_file_name = os.path.join(self.tmpdir, str(time.time()))
        self.submission_builder.cfg['entrypoint_file_name'] = \
                self.entrypoint_file_name
        self.submission_builder.write_json_params = MagicMock()
        self.submission_builder.generate_entrypoint_content = MagicMock(
            return_value='some_content')

    def test_writes_entrypoint(self):
        self.submission_builder.write_entrypoint()
        self.assertEqual(
            self.submission_builder.generate_entrypoint_content.call_args,
            call(params=self.submission_builder.write_json_params\
                 .return_value)
        )
        actual_entrypoint_content = open(self.entrypoint_file_name).read()
        expected_entrypoint_content = \
                self.submission_builder.generate_entrypoint_content.return_value
        self.assertEqual(actual_entrypoint_content, expected_entrypoint_content)

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

class GenerateEntrypointTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_engine_cfg = defaultdict(MagicMock)
        self.submission_builder.cfg['JOB_ENGINE'] = self.job_engine_cfg
        self.params = {('param_%s' % i): ('value_%s' % i) for i in range(3)}

    def test_generates_expected_entrypoint(self):
        expected_entrypoint_content = textwrap.dedent(
            """
            #!/bin/bash
            {job_engine_preamble}
            python -m {job_engine_module} {job_engine_command} {params}
            """
        ).lstrip().format(
            job_engine_preamble=self.job_engine_cfg['entrypoint_preamble'],
            job_engine_module=self.job_engine_cfg['engine_module'],
            job_engine_command='run_job',
            params=self.submission_builder.params_to_cli_args(
                params=self.params)
        )
        self.assertEqual(
            self.submission_builder.generate_entrypoint_content(
                params=self.params),
            expected_entrypoint_content)
