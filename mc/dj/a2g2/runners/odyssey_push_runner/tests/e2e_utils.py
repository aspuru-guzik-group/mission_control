import json
from unittest.mock import patch

from django.conf.urls import url, include
from django.test import TestCase

from ..commands import base as base_command

BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s/' % BASE_PATH, include('missions.urls')),
    url(r'^%s/' % BASE_PATH, include('jobs.urls')),
]

class BaseTestCase(TestCase):
    def setUp(self):
        self.configure_request_client()
        self.request_client_patcher = patch.object(
            base_command, 'request_client', new=self.client)
        self.request_client_patcher.start()

    def configure_request_client(self):
        for method_name in ['patch', 'post']:
            self.patch_client_method_to_use_json(method_name=method_name)

    def patch_client_method_to_use_json(self, method_name=None):
        orig_method = getattr(self.client, method_name)
        def patched_method(path, data=None, **kwargs):
            return orig_method(path, json.dumps(data), 
                               content_type='application/json', **kwargs)
        setattr(self.client, method_name, patched_method)

    def tearDown(self):
        self.request_client_patcher.stop()
