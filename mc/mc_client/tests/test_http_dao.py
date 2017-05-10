import json
import unittest
from unittest.mock import call, MagicMock

from .. import http_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.requests = MagicMock()
        for method in ['post', 'get', 'patch']:
            mock_method = getattr(self.requests, method)
            mock_method.return_value.status_code = 200
        self.dao = http_dao.HttpDao(base_url=self.base_url,
                                    requests=self.requests)

class CreateFlowTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        flow = MagicMock()
        result = self.dao.create_flow(flow=flow)
        self.assertEqual(self.requests.post.call_args,
                         call(self.base_url + 'flows/', json=flow))
        self.assertEqual(result,
                         self.requests.post.return_value.json.return_value)

class GetFlowsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.flow_query_to_query_params = MagicMock()

    def _get_flows(self):
        return self.dao.get_flows(query=self.query)

    def test_generates_query_params(self):
        self._get_flows()
        self.assertEqual(self.dao.flow_query_to_query_params.call_args,
                         call(query=self.query))

    def test_dispatches_to_requests(self):
        result = self._get_flows()
        self.assertEqual(self.requests.get.call_args,
                         call(self.base_url + 'flows/',
                              self.dao.flow_query_to_query_params.return_value))
        self.assertEqual(result,
                         self.requests.get.return_value.json.return_value)

class PatchFlowTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        key = MagicMock()
        patches = MagicMock()
        result = self.dao.patch_flow(key=key, patches=patches)
        self.assertEqual(
            self.requests.patch.call_args,
            call(self.base_url + 'flows/' + key + '/', json=patches)
        )
        self.assertEqual(result,
                         self.requests.patch.return_value.json.return_value)

class CreateJobTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        job = MagicMock()
        result = self.dao.create_job(job=job)
        self.assertEqual(self.requests.post.call_args,
                         call(self.base_url + 'jobs/', json=job))
        self.assertEqual(result,
                         self.requests.post.return_value.json.return_value)

class GetJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.job_query_to_query_params = MagicMock()

    def _get_jobs(self):
        return self.dao.get_jobs(query=self.query)

    def test_generates_query_params(self):
        self._get_jobs()
        self.assertEqual(self.dao.job_query_to_query_params.call_args,
                         call(query=self.query))

    def test_dispatches_to_requests(self):
        result = self._get_jobs()
        self.assertEqual(self.requests.get.call_args,
                         call(self.base_url + 'jobs/',
                              self.dao.job_query_to_query_params.return_value))
        self.assertEqual(result,
                         self.requests.get.return_value.json.return_value)

class PatchJobTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        key = MagicMock()
        patches = MagicMock()
        result = self.dao.patch_job(key=key, patches=patches)
        self.assertEqual(
            self.requests.patch.call_args,
            call(self.base_url + 'jobs/' + key + '/', json=patches)
        )
        self.assertEqual(result,
                         self.requests.patch.return_value.json.return_value)

class FlushTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        result = self.dao.flush_mc_db()
        self.assertEqual(self.requests.get.call_args,
                         call(self.base_url + 'flush/'))
        self.assertEqual(result,
                         self.requests.get.return_value.json.return_value)

class CreateQueueTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        queue = MagicMock()
        result = self.dao.create_queue(queue=queue)
        self.assertEqual(self.requests.post.call_args,
                         call(self.base_url + 'queues/', json=queue))
        self.assertEqual(result,
                         self.requests.post.return_value.json.return_value)

class ClaimQueueItemsTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        queue_key = 'some_key'
        params = {'some': 'param'}
        result = self.dao.claim_queue_items(queue_key=queue_key, params=params)
        self.assertEqual(
            self.requests.post.call_args,
            call(self.base_url + 'queues/' + queue_key + '/claim_items/',
                 data=json.dumps(params))
        )
        self.assertEqual(result,
                         self.requests.post.return_value.json.return_value)

if __name__ == '__main__':
    unittest.main()


