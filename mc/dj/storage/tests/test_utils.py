import os
import tempfile
from django.conf import settings
from django.test import TestCase, override_settings

from .. import utils as storage_utils


class GetStorageBackendTestCase(TestCase):
    def test_gets_storage_backend(self):
        backend = storage_utils.get_storage_backend()
        self.assertTrue(isinstance(backend, storage_utils.FileSystemBackend))

class FileSystemBackendTestCase(TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.backend = storage_utils.FileSystemBackend(base_dir=self.tmpdir)

@override_settings()
class GetBaseDirFromSettingsTestCase(FileSystemBackendTestCase):
    def setUp(self):
        super().setUp()
        attr_name = 'STORAGE_FILESYSTEM_BACKEND_BASEDIR'
        self.mock_base_dir = 'mock_base_dir'
        setattr(settings, attr_name, self.mock_base_dir)

    def test_gets_basedir_from_settings(self):
        base_dir_from_settings = self.backend.get_base_dir_from_settings()
        self.assertEqual(base_dir_from_settings, self.mock_base_dir)

class FileSystemBackendPostDataTestCase(FileSystemBackendTestCase):
    def test_creates_bytes_at_random_uuid(self):
        byte_data = b'some bytes'
        post_result = self.backend.post_data(data=byte_data)
        expected_file_path = os.path.join(self.backend.base_dir,
                                          post_result['key'])
        file_contents = open(expected_file_path, 'rb').read()
        self.assertEqual(file_contents, byte_data)

    def test_encodes_string_as_bytes(self):
        str_data = 'some str'
        post_result = self.backend.post_data(data=str_data)
        expected_file_path = os.path.join(self.backend.base_dir,
                                          post_result['key'])
        file_contents = open(expected_file_path, 'rb').read()
        self.assertEqual(file_contents, str_data.encode('utf-8'))

class FileSystemBackendGetDataTestCase(FileSystemBackendTestCase):
    def test_returns_bytes_at_uuid(self):
        byte_data = b'some bytes'
        post_result = self.backend.post_data(data=byte_data)
        get_result = self.backend.get_data(params={'key': post_result['key']})
        self.assertEqual(get_result, {'data': byte_data})
