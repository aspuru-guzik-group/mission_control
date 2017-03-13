import unittest
from unittest.mock import call, patch

from ...utils import test_utils as job_module_test_utils
from .. import confgen_job_engine


class BaseTestCase(unittest.TestCase,
                   job_module_test_utils.JobEngineTestCaseMixin):
    def setUp(self):
        unittest.TestCase.setUp(self)
        job_module_test_utils.JobEngineTestCaseMixin.setUp(self)
        self.output_dir = 'some_output_dir'

    def generate_engine(self):
        return confgen_job_engine.generate_job_engine()

class GenerateConformersTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job['job_spec']['command'] = 'generate_conformers'

    def decorate_patchers(self):
        self.patchers['confgen_generator'] = patch.object(
            confgen_job_engine.confgen_generator, 'generate_conformers')

    def test_calls_generate_conformers(self):
        self.execute_job(job=self.job, output_dir=self.output_dir)
        self.assertEqual(
            confgen_job_engine.confgen_generator.generate_conformers.call_args,
            call(**{**self.job['job_spec']['kwargs'],
                    'output_dir': self.output_dir})
        )

class ParseCompletedConfgenDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job['job_spec']['command'] = 'parse_completed_confgen_dir'

    def decorate_patchers(self):
        self.patchers['confgen_parser'] = patch.object(
            confgen_job_engine.confgen_parser, 'parse_completed_confgen_dir')

    def test_calls_parse_completed_confgen_dir(self):
        self.execute_job(job=self.job, output_dir=self.output_dir)
        self.assertEqual(
            confgen_job_engine.confgen_parser.\
            parse_completed_confgen_dir.call_args,
            call(
                completed_confgen_dir=self.job['data']['input']['dir_to_parse'],
                output_dir=self.output_dir
            )
        )
