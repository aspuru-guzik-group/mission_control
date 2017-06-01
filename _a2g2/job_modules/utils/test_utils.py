import collections
from unittest.mock import MagicMock


class JobEngineTestCaseMixin(object):
    def setUp(self):
        self.engine = self.generate_engine()
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = self.start_patchers()
        self.job = self.generate_job()

    def generate_engine(self): raise NotImplementedError

    def generate_job(self):
        job = collections.defaultdict(MagicMock)
        job['job_spec'] = collections.defaultdict(MagicMock)
        return job

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
