import unittest
from unittest.mock import call, DEFAULT, Mock, patch
from .storage_client import MissionControlStorageClient


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.storage_client = MissionControlStorageClient(
            base_url=self.base_url)
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def decorate_patchers(self):
        self.patchers['requests'] =  patch.multiple(
            'requests', get=DEFAULT, post=DEFAULT, patch=DEFAULT)

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

class PostTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mocks['requests']['post'].return_value.status_code = 200
        self.data = Mock()
        self.storage_params = Mock()

    def test_makes_post(self):
        self.do_post()
        expected_url = self.base_url + 'post/'
        self.assertEqual(self.mocks['requests']['post'].call_args,
                         call(expected_url, data={
                             'data': self.data,
                             'params': self.storage_params,
                         }))

    def do_post(self):
        return self.storage_client.post_data(data=self.data,
                                             storage_params=self.storage_params)

    def test_returns_response_data(self):
        response = {'key': Mock()}
        self.mocks['requests']['post'].return_value.json.return_value = response
        result = self.do_post()
        self.assertEqual(result, response)

class GetTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mocks['requests']['get'].return_value.status_code = 200
        self.storage_params = {'key': Mock()}

    def test_makes_get(self):
        self.do_get()
        expected_url = self.base_url + 'get/'
        self.assertEqual(
            self.mocks['requests']['get'].call_args,
            call(expected_url, data={'params': self.storage_params}))

    def do_get(self):
        return self.storage_client.get_data(storage_params=self.storage_params)

    def test_returns_response_data(self):
        response = {'data': Mock()}
        self.mocks['requests']['get'].return_value.json.return_value = response
        result = self.do_get()
        self.assertEqual(result, response)

if __name__ == '__main__':
    unittest.main()


