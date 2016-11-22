import unittest
from unittest.mock import call, DEFAULT, patch
from jobs.models import Job
from .job_spec_client import MissionControlJobSpecClient


class JobSpecClientBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.client = MissionControlJobSpecClient(
            base_url=self.base_url
        )
        self.patchers = {
            'requests': patch.multiple('requests', get=DEFAULT, post=DEFAULT,
                                       patch=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

class FetchJobSpecsTestCase(JobSpecClientBaseTestCase):
    def test_makes_get(self):
        self.client.fetch_job_specs()
        expected_url = self.base_url + 'jobs/'
        expected_get_call = call(expected_url)
        self.assertEqual(self.mocks['requests']['get'].call_args,
                         expected_get_call)

    def test_filters(self):
        filters = {'filter1': 'junk'}
        self.client.fetch_job_specs(filters=filters)
        expected_url = self.base_url + 'jobs/'
        self.assertEqual(self.mocks['requests']['get'].call_args,
                         call(expected_url, filters))

    def test_handles_populated_response(self):
        jobs = [{'uuid': i} for i in range(3)]
        self.mocks['requests']['get'].return_value.json.return_value = jobs
        fetched_job_specs = self.client.fetch_job_specs()
        self.assertEqual(fetched_job_specs, jobs)

    def test_handles_empty_response(self):
        jobs = []
        self.mocks['requests']['get'].return_value.json.return_value = jobs
        fetched_job_specs = self.client.fetch_job_specs()
        self.assertEqual(fetched_job_specs, jobs)

class FetchClaimableJobSpecsTest(JobSpecClientBaseTestCase):
    def test_calls_fetch_job_specs_with_filters(self):
        with patch.multiple(self.client, fetch_job_specs=DEFAULT):
            self.client.fetch_claimable_job_specs()
            self.assertEqual(self.client.fetch_job_specs.call_args,
                             call(filters={'status': Job.STATUSES.PENDING}))

class ClaimJobSpecsTestCase(JobSpecClientBaseTestCase):
    def setUp(self):
        super().setUp()
        self.uuids = [i for i in range(3)]

    def test_makes_post(self):
        self.client.claim_job_specs(uuids=self.uuids)
        expected_url = self.base_url + 'claim_jobs/'
        expected_post_call = call(expected_url, {'uuids': self.uuids})
        self.assertEqual(self.mocks['requests']['post'].call_args,
                         expected_post_call)

    def test_handles_populated_response(self):
        self.mocks['requests']['post'].return_value.json.return_value = {
            _uuid: True for _uuid in self.uuids
        }
        result = self.client.claim_job_specs(uuids=self.uuids)
        expected_result = {_uuid: True for _uuid in self.uuids}
        self.assertEqual(result, expected_result)

    def test_handles_empty_response(self):
        self.mocks['requests']['post'].return_value.json.return_value = {}
        result = self.client.claim_job_specs(uuids=self.uuids)
        expected_result = {}
        self.assertEqual(result, expected_result)

class UpdateJobSpecsTestCase(JobSpecClientBaseTestCase):
    def setUp(self):
        super().setUp()
        self.updates_by_uuid = {
            'uuid_%s' % i: {'key1': 'value1', 'key2': 'value2'}
            for i in range(3)
        }

    def test_makes_patch_calls(self):
        self.client.update_job_specs(updates_by_uuid=self.updates_by_uuid)
        expected_patch_calls = [
            call(self.base_url + 'jobs/' + _uuid + '/', updates_for_uuid)
            for _uuid, updates_for_uuid in self.updates_by_uuid.items()
        ]
        self.assertEqual(
            sorted(self.mocks['requests']['patch'].call_args_list,
                   key=lambda _call_as_pair: _call_as_pair[0][0]),
            sorted(expected_patch_calls,
                   key=lambda _call_as_triple: _call_as_triple[1][0]))

    def test_returns_combined_patch_results(self):
        self.mocks['requests']['patch'].return_value.json.return_value = {}
        result = self.client.update_job_specs(
            updates_by_uuid=self.updates_by_uuid)
        expected_result = {_uuid: {} for _uuid in  self.updates_by_uuid}
        self.assertEqual(result, expected_result)

if __name__ == '__main__':
    unittest.main()


