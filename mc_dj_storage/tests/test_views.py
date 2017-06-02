import io
import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings
from unittest.mock import call, patch

from .. import views as _views
from .. import urls as _urls


BASE_PATH = 'storage'
BASE_URL = '/%s' % BASE_PATH
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_urls.__name__, namespace='storage')),
]

class _BaseTestCase(TestCase):
    def setUp(self):
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.data = b'some data'
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
        self.expected_backend.post_data.return_value = 'some backend response'
        self.result = self.do_post()

    def do_post(self, data_to_post=None):
        post_data = data_to_post or {
            'data': io.BytesIO(self.data),
            'params': io.StringIO(json.dumps(self.params))
        }
        response = self.client.post(self.url, data=post_data)
        if not str(response.status_code).startswith('2'):
            raise Exception("Bad response: %s" % response)
        return response.json()

    def test_gets_backend_for_params(self):
        self.assertEqual(_views.storage_utils.get_storage_backend.call_args,
                         call(params=self.params))

    def test_calls_backend_post_data_with_data(self):
        call_args = self.expected_backend.post_data.call_args
        data_arg = call_args[1]['data']
        params_arg = call_args[1]['params']
        self.assertEqual(data_arg.name, 'data')
        self.assertEqual(params_arg, self.params)

    def test_returns_result(self):
        expected_params = self.expected_backend.post_data.return_value
        self.assertEqual(
            self.result,
            {
                'params': expected_params,
                'url': _views.get_url_for_params(params=expected_params)
            }
        )

class TestGetDataTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.url = BASE_URL + 'get/'
        self.expected_backend.get_data.return_value = b'some backend response'
        self.response = self.do_get()

    def do_get(self):
        response = self.client.get(self.url, data={
            'params': json.dumps(self.params)})
        if not str(response.status_code).startswith('2'):
            raise Exception("Bad response: %s" % response)
        return response

    def test_gets_backend_for_params(self):
        self.assertEqual(_views.storage_utils.get_storage_backend.call_args,
                         call(params=self.params))

    def test_calls_backend_get_data(self):
        self.assertEqual(self.expected_backend.get_data.call_args,
                         call(params=self.params))

    def test_returns_raw_data(self):
        self.assertEqual(self.response.content,
                         self.expected_backend.get_data.return_value)
