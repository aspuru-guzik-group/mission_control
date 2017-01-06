import unittest
from unittest.mock import call, DEFAULT, patch
from .flow_client import MissionControlFlowClient


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.client = MissionControlFlowClient(base_url=self.base_url)
        self.patchers = {
            'requests': patch.multiple('requests', get=DEFAULT, post=DEFAULT,
                                       patch=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

class FetchFlowsTestCase(BaseTestCase):
    def test_makes_get(self):
        self.client.fetch_flows()
        expected_url = self.base_url + 'flows/'
        expected_get_call = call(expected_url)
        self.assertEqual(self.mocks['requests']['get'].call_args,
                         expected_get_call)

    def test_query_params(self):
        query_params = {'filter1': 'junk'}
        self.client.fetch_flows(query_params=query_params)
        expected_url = self.base_url + 'flows/'
        self.assertEqual(self.mocks['requests']['get'].call_args,
                         call(expected_url, query_params))

    def test_handles_populated_response(self):
        flows = [{'uuid': i} for i in range(3)]
        self.mocks['requests']['get'].return_value.json.return_value = flows
        fetched_flows = self.client.fetch_flows()
        self.assertEqual(fetched_flows, flows)

    def test_handles_empty_response(self):
        flows = []
        self.mocks['requests']['get'].return_value.json.return_value = flows
        fetched_flows = self.client.fetch_flows()
        self.assertEqual(fetched_flows, flows)

class ClaimFlowsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.uuids = [i for i in range(3)]

    def test_makes_post(self):
        self.client.claim_flows(uuids=self.uuids)
        expected_url = self.base_url + 'claim_flows/'
        expected_post_call = call(expected_url, {'uuids': self.uuids})
        self.assertEqual(self.mocks['requests']['post'].call_args,
                         expected_post_call)

    def test_handles_populated_response(self):
        self.mocks['requests']['post'].return_value.json.return_value = {
            _uuid: True for _uuid in self.uuids
        }
        result = self.client.claim_flows(uuids=self.uuids)
        expected_result = {_uuid: True for _uuid in self.uuids}
        self.assertEqual(result, expected_result)

    def test_handles_empty_response(self):
        self.mocks['requests']['post'].return_value.json.return_value = {}
        result = self.client.claim_flows(uuids=self.uuids)
        expected_result = {}
        self.assertEqual(result, expected_result)

class UpdateFlowsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.updates_by_uuid = {
            'uuid_%s' % i: {'key1': 'value1', 'key2': 'value2'}
            for i in range(3)
        }

    def test_makes_patch_calls(self):
        self.client.update_flows(updates_by_uuid=self.updates_by_uuid)
        expected_patch_calls = [
            call(self.base_url + 'flows/' + _uuid + '/', updates_for_uuid)
            for _uuid, updates_for_uuid in self.updates_by_uuid.items()
        ]
        self.assertEqual(
            sorted(self.mocks['requests']['patch'].call_args_list,
                   key=lambda _call_as_pair: _call_as_pair[0][0]),
            sorted(expected_patch_calls,
                   key=lambda _call_as_triple: _call_as_triple[1][0]))

    def test_returns_combined_patch_results(self):
        self.mocks['requests']['patch'].return_value.json.return_value = {}
        result = self.client.update_flows(
            updates_by_uuid=self.updates_by_uuid)
        expected_result = {_uuid: {} for _uuid in  self.updates_by_uuid}
        self.assertEqual(result, expected_result)

class CreateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = {'data': 'some data'}

    def test_makes_post_call(self):
        self.client.create_flow(flow=self.flow)
        self.assertEqual(self.mocks['requests']['post'].call_args,
                         call(self.base_url + 'flows/', data=self.flow))

    def test_returns_post_result(self):
        mock_result = {'some': 'result'}
        self.mocks['requests']['post'].return_value.json.return_value = \
                mock_result
        result = self.client.create_flow(flow=self.flow)
        self.assertEqual(result, mock_result)

if __name__ == '__main__':
    unittest.main()


