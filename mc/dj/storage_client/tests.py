import unittest
from unittest.mock import call, Mock, patch
from . import storage_client


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.storage_client = storage_client.MissionControlStorageClient(
            base_url=self.base_url)

    def decorate_patchers(self):
        self.patchers['requests'] = patch.object(storage_client, 'requests')

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

class PostTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mocks['requests'].post.return_value.status_code = 200
        self.data = Mock()
        self.storage_params = Mock()
        self.storage_client.serialize_storage_params = Mock()

    def test_makes_post(self):
        self.do_post()
        expected_url = self.base_url + 'post/'
        self.assertEqual(
            self.mocks['requests'].post.call_args,
            call(expected_url, files=[
                ('data', ('data', self.data, 'application/octet-stream')),
                ('params', ('params',
                            self.storage_client.serialize_storage_params\
                                .return_value,
                            'application/json'))
            ])
        )

    def do_post(self):
        return self.storage_client.post_data(data=self.data,
                                             storage_params=self.storage_params)

    def test_returns_response_params(self):
        response = {'params': {'key': Mock()}}
        self.mocks['requests'].post.return_value.json.return_value = response
        result = self.do_post()
        self.assertEqual(result, response['params'])

class GetTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mocks['requests'].get.return_value.status_code = 200
        self.storage_params = {'key': Mock()}
        self.storage_client.serialize_storage_params = Mock()

    def test_makes_get(self):
        self.do_get()
        expected_url = self.base_url + 'get/'
        self.assertEqual(
            self.mocks['requests'].get.call_args,
            call(
                expected_url,
                params={
                    'params': self.storage_client.serialize_storage_params\
                        .return_value
                },
                stream=True
            )
        )

    def do_get(self):
        return self.storage_client.get_data(storage_params=self.storage_params)

    def test_returns_raw_data(self):
        result = self.do_get()
        self.assertEqual(result,
                         self.mocks['requests'].get.return_value.content)

class SerializeStorageParamsTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['json'] = patch.object(storage_client, 'json')

    def test_serialize_storage_params(self):
        storage_params = Mock()
        result = self.storage_client.serialize_storage_params(storage_params)
        self.assertEqual(self.mocks['json'].dumps.call_args,
                         call(storage_params))
        self.assertEqual(result, self.mocks['json'].dumps.return_value)

if __name__ == '__main__':
    unittest.main()


