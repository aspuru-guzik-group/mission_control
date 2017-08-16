import unittest

from mc.houston.tests import utils as _houston_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = _houston_test_utils.generate_test_houston()
        self.jobs = self.create_jobs()

    def _generate_cfg(self):
        return _houston_test_utils.generate_test_cfg()

    def create_jobs(self):
        jobs = [
            self.Job(
                job_type=('job_type_%s' % i),
                job_params={'param': ('param_%s' % i)}
            )
            for i in range(3)
        ]
        self.session.add_all(jobs)
        self.session.commit()
        return jobs

    def _delete_items(self, **kwargs):
        return self.houston.run_command('delete_items', **kwargs)

    @property
    def Job(self): return self.houston.utils.db.models.Job

    @property
    def session(self): return self.houston.utils.db.session


class QueryItemsTestCase(BaseTestCase):
    def test_returns_query_results(self):
        query = {
            'filters': [
                {'field': 'job_type', 'op': 'IN',
                 'arg': ['job_type_0', 'job_type_2']}
            ]
        }
        item_type = 'Job'
        results = self._delete_items(item_type=item_type, query=query)
        expected_results = {'num_deleted': 2}
        self.assertEqual(results, expected_results)
        jobs_in_db = self.session.query(self.Job).all()
        self.assertEqual(len(jobs_in_db), 1)
        self.assertEqual(jobs_in_db[0].job_type, 'job_type_1')
