import unittest

from mc.houston.tests import utils as _houston_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = _houston_test_utils.generate_test_houston()
        self.job_type = 'some_job_type'
        self.job_params = {'params': 'some_job_params'}

    def _generate_cfg(self):
        return _houston_test_utils.generate_test_cfg()

    def _get_or_create_job(self, flow_spec=None):
        flow_spec = flow_spec or {}
        return self.houston.run_command(
            'get_or_create_job',
            job_type=self.job_type,
            job_params=self.job_params
        )

    @property
    def session(self): return self.houston.utils.db.session

    @property
    def Job(self): return self.houston.utils.db.models.Job


class GetOrCreateJobTestCase(BaseTestCase):
    def test_creates_job_if_not_exists(self):
        job_dict = self._get_or_create_job()
        jobs_in_db = self.session.query(self.Job).all()
        self.assertEqual(len(jobs_in_db), 1)
        self.assertEqual(jobs_in_db[0].to_dict(), job_dict)

    def test_returns_existing_job(self):
        existing_job = self.Job(
            job_type=self.job_type, job_params=self.job_params)
        with self.session.begin(subtransactions=True):
            self.session.add(existing_job)
        jobs_in_db = self.session.query(self.Job).all()
        self.assertEqual(len(jobs_in_db), 1)
        job_dict = self._get_or_create_job()
        jobs_in_db = self.session.query(self.Job).all()
        self.assertEqual(len(jobs_in_db), 1)
        self.assertEqual(existing_job.to_dict(), job_dict)
