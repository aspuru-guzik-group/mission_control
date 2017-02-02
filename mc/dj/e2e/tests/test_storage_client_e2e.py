import tempfile
from django.conf.urls import url, include
from django.conf import settings
from django.test import TestCase, override_settings

from storage_client.storage_client import MissionControlStorageClient

from storage import urls as _storage_urls

BASE_PATH = 'storage/'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_storage_urls.__name__)),
]

@override_settings(ROOT_URLCONF=__name__)
class StorageClientE2ETestCase(TestCase):
    def setUp(self):
        self.setup_storage_basedir()
        self.storage_client = MissionControlStorageClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )
        self.data = 'some data'
        self.storage_params = 'some params'

    def setup_storage_basedir(self):
        self.storage_base_dir = tempfile.mkdtemp()
        settings_attr = 'STORAGE_FILESYSTEM_BACKEND_BASEDIR'
        setattr(settings, settings_attr, self.storage_base_dir)

    def test_storage(self):
        updated_params = self.storage_client.post_data(
            data=self.data, storage_params=self.storage_params)
        get_result = self.storage_client.get_data(storage_params=updated_params)
        self.assertEqual(get_result['data'], self.data)
