import unittest

from mc.houston.houston import Houston
from mc.houston.tests import utils as _houston_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = Houston(cfg=self._generate_cfg())
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

    def _query_items(self, **kwargs):
        return self.houston.run_command('query_items', **kwargs)

    @property
    def Job(self): return self.houston.utils.db.models.Job

    @property
    def session(self): return self.houston.utils.db.session


class QueryItemsTestCase(BaseTestCase):
    def test_returns_query_results(self):
        query = {
            'filters': [
                {'field': 'job_type', 'op': '=', 'arg': 'job_type_1'}
            ]
        }
        item_type = 'Job'
        results = self._query_items(item_type=item_type, query=query)
        expected_results = [self.jobs[1].to_dict()]
        self.assertEqual(results, expected_results)
