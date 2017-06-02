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
                                    requests=self.requests,
                                    urls=MagicMock())
        self.item_type = MagicMock()

class CreateItemTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        kwargs = MagicMock()
        result = self.dao.create_item(item_type=self.item_type, kwargs=kwargs)
        self.assertEqual(self.requests.post.call_args,
                         call(self.dao.urls[self.item_type], json=kwargs))
        self.assertEqual(result,
                         self.requests.post.return_value.json.return_value)

class GetItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.query_to_query_params = MagicMock()

    def _get_items(self):
        return self.dao.get_items(item_type=self.item_type, query=self.query)

    def test_generates_query_params(self):
        self._get_items()
        self.assertEqual(self.dao.query_to_query_params.call_args,
                         call(query=self.query))

    def test_dispatches_to_requests(self):
        result = self._get_items()
        self.assertEqual(self.requests.get.call_args,
                         call(self.dao.urls[self.item_type],
                              self.dao.query_to_query_params.return_value))
        self.assertEqual(result,
                         self.requests.get.return_value.json.return_value)

class PatchItemTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        key = MagicMock()
        patches = MagicMock()
        result = self.dao.patch_item(item_type=self.item_type, key=key,
                                     patches=patches)
        self.assertEqual(
            self.requests.patch.call_args,
            call(self.dao.urls[self.item_type] + key + '/', json=patches)
        )
        self.assertEqual(result,
                         self.requests.patch.return_value.json.return_value)

class FlushTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        result = self.dao.flush_mc_db()
        self.assertEqual(self.requests.get.call_args,
                         call(self.dao.urls['flush']))
        self.assertEqual(result,
                         self.requests.get.return_value.json.return_value)

class ClaimQueueItemsTestCase(BaseTestCase):
    def test_dispatches_to_requests(self):
        queue_key = 'some_key'
        result = self.dao.claim_queue_items(queue_key=queue_key)
        expected_url = '{queues_root}{queue_key}/claim_items/'.format(
            queues_root=self.dao.urls['Queue'], queue_key=queue_key)
        self.assertEqual(self.requests.post.call_args, call(expected_url))
        self.assertEqual(result,
                         self.requests.post.return_value.json.return_value)

if __name__ == '__main__':
    unittest.main()


