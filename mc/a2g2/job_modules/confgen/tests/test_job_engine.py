import collections
import unittest
from unittest.mock import call, patch, MagicMock

from .. import job_engine as confgen_job_engine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = confgen_job_engine.ConfgenJobEngine()
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = self.start_patchers()
        self.job = self.generate_job()

    def generate_job(self):
        job = collections.defaultdict(MagicMock)
        job['job_spec'] = collections.defaultdict(MagicMock)

    def decorate_patchers(self): pass

    def start_patchers(self):
        mocks = {key: patcher.start() for key, patcher in self.patchers.items()}
        return mocks

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers()

    def stop_patchers(self):
        for patcher in self.patchers.values(): patcher.stop()

    def execute_job(self, job=None, **kwargs):
        self.engine.execute_job(job=job, **kwargs)

class GenerateConformersTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job['command'] = 'generate_conformers'
        self.output_dir = 'some_output_dir'

    def decorate_patchers(self):
        self.patchers['conformer_generator'] = patch.object(
            confgen_job_engine.conformer_generator, 'generate_conformers')

    def test_calls_generate_conformers(self):
        self.execute_job(job=self.job, output_dir=self.output_dir)
        self.assertEqual(
            confgen_job_engine.conformer_generator.generate_conformers.call_args,
            call(**{**self.job['job_spec']['kwargs'],
                    'output_dir': self.output_dir})
        )




