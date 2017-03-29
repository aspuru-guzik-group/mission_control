from collections import defaultdict
import json
import os
import tempfile
import unittest
from unittest.mock import call, MagicMock, patch

from .. import submission_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.submission = self.generate_submission()
        self.submission_runner = submission_runner.SubmissionRunner(
            submission=self.submission)

    def generate_submission(self):
        self.job = self.generate_job()
        self.cfg = self.generate_cfg()
        submission_dir = tempfile.mkdtemp()
        self.write_json_files(submission_dir=submission_dir)
        submission = {
            'dir': submission_dir,
            'outputs_dir': os.path.join(submission_dir, 'outputs'),
        }
        return submission

    def generate_job(self):
        return defaultdict(MagicMock)

    def generate_cfg(self):
        return defaultdict(MagicMock)

    def write_json_files(self, submission_dir=None):
        for param in ['job', 'cfg']:
            with open(os.path.join(submission_dir, param + '.json'), 'w') as f:
                json.dump(getattr(self, param), f)

class RunSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.expected_output_dir = os.path.join(
            self.submission['outputs_dir'], 'confgen_outputs')
        self.submission_runner.get_confgen_params = MagicMock()
        self.submission_runner.generate_conformers = MagicMock()
        self.submission_runner.run_submission()

    def test_creates_output_dir(self):
        self.assertTrue(os.path.exists(self.expected_output_dir))

    def test_calls_generate_conformers(self):
        self.assertEqual(
            self.submission_runner.generate_conformers.call_args,
            call(confgen_params=\
                 self.submission_runner.get_confgen_params.return_value,
                 output_dir=self.expected_output_dir
                )
        )

class GetConfgenParamsTestCase(BaseTestCase):
    def generate_job(self):
        job = {
            'job_spec': {
                'job_params': {
                    'confgen_params': {
                        'param': 'value'
                    }
                }
            }
        }
        return job

    def test_gets_confgen_params_from_job(self):
        confgen_params = self.submission_runner.get_confgen_params()
        self.assertEqual(confgen_params,
                         self.job['job_spec']['job_params']['confgen_params'])

class GenerateConformersTestCase(BaseTestCase):
    @patch.object(submission_runner, 'ConformerGenerator')
    def test_dispatches_to_conformer_generator(self, MockConformerGenerator):
        confgen_params = defaultdict(MagicMock)
        self.submission_runner.generate_conformers(
            confgen_params=confgen_params)
        self.assertEqual(MockConformerGenerator.call_args,
                         call(**confgen_params))
        self.assertEqual(
            MockConformerGenerator.return_value.generate_conformers.call_args,
            call())
