import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings
from unittest.mock import call, patch

from .. import views as _views
from .. import urls as _urls

BASE_PATH = 'storage'
BASE_URL = '/%s' % BASE_PATH
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_urls.__name__)),
]

class _BaseTestCase(TestCase):
    def setUp(self):
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.data = 'some data'
        self.params = 'some params'
        self.expected_backend = _views.storage_utils.get_storage_backend\
                .return_value

    def decorate_patchers(self):
        self.patchers['storage_utils'] = patch.object(_views, 'storage_utils')

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

BaseTestCase = override_settings(ROOT_URLCONF=__name__)(_BaseTestCase)


class TestPostDataTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.url = BASE_URL + 'post/'
        self.expected_backend.post.return_value = 'some backend response'
        self.result = self.do_post()

    def do_post(self, data_to_post=None):
        post_data = data_to_post or {'data': self.data,
                                     'params': json.dumps(self.params)}
        response = self.client.post(self.url, data=post_data)
        if not str(response.status_code).startswith('2'):
            raise Exception("Bad response: %s" % response)
        return response.json()

    def test_gets_backend_for_params(self):
        self.assertEqual(_views.storage_utils.get_storage_backend.call_args,
                         call(params=self.params))

    def test_calls_backend_post(self):
        self.assertEqual(self.expected_backend.post.call_args,
                         call(data=self.data, params=self.params))

    def test_returns_result(self):
        self.assertEqual(
            self.result,
            {'storage_uri': self.expected_backend.post.return_value}
        )


class TestGetDataTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.url = BASE_URL + 'get/'
        self.expected_backend.get.return_value = 'some backend response'
        self.result = self.do_get()

    def do_get(self):
        response = self.client.get(self.url, data={
            'params': json.dumps(self.params)})
        if not str(response.status_code).startswith('2'):
            raise Exception("Bad response: %s" % response)
        return response.json()

    def test_gets_backend_for_params(self):
        self.assertEqual(_views.storage_utils.get_storage_backend.call_args,
                         call(params=self.params))

    def test_calls_backend_get(self):
        self.assertEqual(self.expected_backend.get.call_args,
                         call(params=self.params))

    def test_returns_result(self):
        self.assertEqual(
            self.result,
            {'data': self.expected_backend.get.return_value}
        )
