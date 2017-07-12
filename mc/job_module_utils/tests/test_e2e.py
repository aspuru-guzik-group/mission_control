import tempfile
import unittest

from ..job_module_command_dispatcher import JobModuleCommandDispatcher
from ..jobdir_builder import JobdirBuilder
from ..jobdir_runner import JobdirRunner

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.db = {}
        self.job_module = self._generate_job_module(db=self.db)
        self.command_dispatcher = JobModuleCommandDispatcher(
            load_job_module_fn=self._load_job_module)
        self.job = {'some': 'job'}
        self.cfg = {'some': 'cfg'}

    def _generate_job_module(self, db=None):
        class MyJobModule(object):
            def build_jobdir(self_, job=None, cfg=None, output_dir=None):
                builder = JobdirBuilder(job=job, cfg=cfg, output_dir=output_dir)
                return builder.build_jobdir()

            def run_jobdir(self_, job_spec=None, job=None, cfg=None):
                runner = JobdirRunner(job_spec=job_spec, job=job, cfg=cfg,
                                      run_jobdir_fn=self_._run_jobdir)
                runner.run_jobdir()

            def _run_jobdir(self_, job_spec=None, job=None, cfg=None,
                            utils=None):
                db['job_spec'] = job_spec
                db['job'] = job
                db['cfg'] = cfg

        return MyJobModule()

    def _load_job_module(self, job=None, cfg=None):
        return self.job_module

class E2ETestCase(BaseTestCase):
    def test_builds_and_runs_dir(self):
        output_dir = tempfile.mkdtemp()
        job_spec = self.command_dispatcher.build_jobdir(
            job=self.job, cfg=self.cfg, output_dir=output_dir)
        self.command_dispatcher.run_jobdir(jobdir=job_spec['dir'], cfg=self.cfg)
        expected_db = {
            'job_spec': job_spec,
            'job': self.job,
            'cfg': self.cfg
        }
        self.assertEqual(self.db, expected_db)
