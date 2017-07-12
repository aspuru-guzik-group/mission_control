import unittest
from unittest.mock import call, MagicMock

from .. import jobdir_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.mod_mocks = {}
        self.job_spec = MagicMock()
        self.job = MagicMock()
        self.cfg = MagicMock()
        self.run_jobdir_fn = MagicMock()
        self.runner = jobdir_runner.JobdirRunner(
            job_spec=self.job_spec,
            job=self.job,
            cfg=self.cfg,
            run_jobdir_fn=self.run_jobdir_fn
        )

class RunJobdirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.run_jobdir()

    def test_dispatches_to_run_jobdir_fn(self):
        self.assertEqual(
            self.runner.run_jobdir_fn.call_args,
            call(job_spec=self.runner.job_spec, job=self.runner.job,
                 cfg=self.runner.cfg, utils=self.runner.utils)
        )
