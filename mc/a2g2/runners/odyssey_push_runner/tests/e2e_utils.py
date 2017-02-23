from unittest.mock import patch

from django.conf.urls import url, include
from django.test import TestCase

from mc.mc_utils import test_utils

from ..commands import base as base_command

BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s/' % BASE_PATH, include('missions.urls')),
]

class BaseTestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client(request_client=self.client)
        self.request_client_patcher = patch.object(
            base_command, 'request_client', new=self.client)
        self.request_client_patcher.start()

    def tearDown(self):
        self.request_client_patcher.stop()
