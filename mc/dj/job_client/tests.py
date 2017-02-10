import json
import unittest
from unittest.mock import call, DEFAULT, Mock, patch
from jobs.models import JobStatuses
from .job_client import MissionControlJobClient


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.job_client = MissionControlJobClient(
            base_url=self.base_url
        )
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def decorate_patchers(self):
        self.patchers['requests'] =  patch.multiple(
            'requests', get=DEFAULT, post=DEFAULT, patch=DEFAULT)

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

class FetchJobsTestCase(BaseTestCase):
    def test_makes_get(self):
        self.job_client.fetch_jobs()
        expected_url = self.base_url + 'jobs/'
        expected_get_call = call(expected_url)
        self.assertEqual(self.mocks['requests']['get'].call_args,
                         expected_get_call)

    def test_query_params(self):
        query_params = {'filter1': 'junk'}
        self.job_client.fetch_jobs(query_params=query_params)
        expected_url = self.base_url + 'jobs/'
        self.assertEqual(self.mocks['requests']['get'].call_args,
                         call(expected_url, query_params))

    def test_handles_populated_response(self):
        jobs = [{'uuid': i} for i in range(3)]
        self.mocks['requests']['get'].return_value.json.return_value = jobs
        fetched_jobs = self.job_client.fetch_jobs()
        self.assertEqual(fetched_jobs, jobs)

    def test_handles_empty_response(self):
        jobs = []
        self.mocks['requests']['get'].return_value.json.return_value = jobs
        fetched_jobs = self.job_client.fetch_jobs()
        self.assertEqual(fetched_jobs, jobs)

class FetchClaimableJobsTest(BaseTestCase):
    def test_calls_fetch_jobs_with_query_params(self):
        with patch.multiple(self.job_client, fetch_jobs=DEFAULT):
            self.job_client.fetch_claimable_jobs()
            self.assertEqual(
                self.job_client.fetch_jobs.call_args,
                call(query_params={'status': JobStatuses.PENDING.name}))

class FetchJobByUUIDTestCase(BaseTestCase):
    def decorate_patchers(self):
        super().decorate_patchers()
        self.patchers['client'] = patch.multiple(self.job_client,
                                                 fetch_jobs=DEFAULT)

    def test_wraps_fetch_jobs(self):
        mock_fetch_jobs = self.mocks['client']['fetch_jobs']
        mock_fetch_jobs.return_value = [Mock()]
        uuid = 'some uuid'
        result = self.job_client.fetch_job_by_uuid(uuid=uuid)
        self.assertEqual(mock_fetch_jobs.call_args,
                         call(query_params={'uuid': uuid}))
        self.assertEqual(result, mock_fetch_jobs.return_value[0])

class ClaimJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.uuids = [i for i in range(3)]

    def test_makes_post(self):
        self.job_client.claim_jobs(uuids=self.uuids)
        expected_url = self.base_url + 'claim_jobs/'
        expected_post_call = call(expected_url, json={'uuids': self.uuids})
        self.assertEqual(self.mocks['requests']['post'].call_args,
                         expected_post_call)

    def test_handles_populated_response(self):
        self.mocks['requests']['post'].return_value.json.return_value = {
            _uuid: True for _uuid in self.uuids
        }
        result = self.job_client.claim_jobs(uuids=self.uuids)
        expected_result = {_uuid: True for _uuid in self.uuids}
        self.assertEqual(result, expected_result)

    def test_handles_empty_response(self):
        self.mocks['requests']['post'].return_value.json.return_value = {}
        result = self.job_client.claim_jobs(uuids=self.uuids)
        expected_result = {}
        self.assertEqual(result, expected_result)

class UpdateJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.updates_by_uuid = {
            'uuid_%s' % i: {'key1': 'value1', 'key2': 'value2'}
            for i in range(3)
        }

    def test_makes_patch_calls(self):
        self.job_client.update_jobs(updates_by_uuid=self.updates_by_uuid)
        expected_patch_calls = [
            call(self.base_url + 'jobs/' + _uuid + '/',
                 json=self.job_client.format_job_kwargs(updates_for_uuid))
            for _uuid, updates_for_uuid in self.updates_by_uuid.items()
        ]
        self.assertEqual(
            sorted(self.mocks['requests']['patch'].call_args_list,
                   key=lambda _call_as_pair: _call_as_pair[0][0]),
            sorted(expected_patch_calls,
                   key=lambda _call_as_triple: _call_as_triple[1][0]))

    def test_returns_combined_patch_results(self):
        self.mocks['requests']['patch'].return_value.json.return_value = {}
        result = self.job_client.update_jobs(
            updates_by_uuid=self.updates_by_uuid)
        expected_result = {_uuid: {} for _uuid in  self.updates_by_uuid}
        self.assertEqual(result, expected_result)

class CreateJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_kwargs = {'data': 'some data'}
        self.mocks['requests']['post'].return_value.status_code = 200

    def test_makes_post_call(self):
        self.job_client.create_job(job_kwargs=self.job_kwargs)
        self.assertEqual(
            self.mocks['requests']['post'].call_args,
            call(self.base_url + 'jobs/', json={
                **self.job_kwargs,
                'data': json.dumps(self.job_kwargs['data']),
            })
        )

    def test_returns_post_result(self):
        mock_result = {'some': 'result'}
        self.mocks['requests']['post'].return_value.json.return_value = \
                mock_result
        result = self.job_client.create_job(job_kwargs=self.job_kwargs)
        self.assertEqual(result, mock_result)

if __name__ == '__main__':
    unittest.main()


