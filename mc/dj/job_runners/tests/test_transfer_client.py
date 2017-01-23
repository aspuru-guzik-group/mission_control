import os
import unittest
from unittest.mock import call, Mock

from ..transfer_client import TransferClient

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.mkdtemp = Mock(return_value='/mock_tmp_dir')
        self.transfer_client = TransferClient(mkdtemp=self.mkdtemp)
        self.expected_staging_dir = os.path.join(
            self.mkdtemp.return_value, self.transfer_client.staging_dir_name)
        self.src = Mock()
        self.dest = Mock()

class CopyTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.transfer_client.download = Mock()
        self.transfer_client.upload = Mock()
        self.result = self.transfer_client.copy(src=self.src, dest=self.dest)

    def test_downloads_src_to_staging_dir(self):
        self.assertEqual(self.transfer_client.download.call_args,
                         call(src=self.src, dest=self.expected_staging_dir))

    def test_uploads_staging_dir_to_dest(self):
        self.assertEqual(self.transfer_client.upload.call_args,
                         call(src=self.expected_staging_dir, dest=self.dest))

    def test_returns_upload_output(self):
        self.assertEqual(self.result, self.transfer_client.upload.return_value)

class DownloadTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.transfer_client.get_subclient_for_transfer_uri = Mock()
        self.result = self.transfer_client.download(src=self.src,
                                                    dest=self.dest)

    def test_gets_subclient_for_src(self):
        self.assertEqual(
            self.transfer_client.get_subclient_for_transfer_uri.call_args,
            call(transfer_uri=self.src)
        )

    def test_returns_subclient_download(self):
        subclient = self.transfer_client.get_subclient_for_transfer_uri\
                .return_value
        self.assertEqual(subclient.download.call_args,
                         call(src=self.src, dest=self.dest))
        self.assertEqual(self.result, subclient.download.return_value)

class UploadTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.transfer_client.get_subclient_for_transfer_uri = Mock()
        self.result = self.transfer_client.upload(src=self.src, dest=self.dest)

    def test_gets_subclient_for_dest(self):
        self.assertEqual(
            self.transfer_client.get_subclient_for_transfer_uri.call_args,
            call(transfer_uri=self.dest)
        )

    def test_returns_subclient_upload(self):
        subclient = self.transfer_client.get_subclient_for_transfer_uri\
                .return_value
        self.assertEqual(subclient.upload.call_args,
                         call(src=self.src, dest=self.dest))
        self.assertEqual(self.result, subclient.upload.return_value)

class GetSubclientForTransferUriTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.storage_key = 'some_storage'
        self.object_key = 'some_object'
        self.transfer_uri = ':'.join([self.storage_key, self.object_key])
        self.subclient = Mock()
        self.transfer_client.register_subclient(subclient=self.subclient,
                                                key=self.storage_key)

    def test_gets_subclient_for_transfer_uri(self):
        result = self.transfer_client.get_subclient_for_transfer_uri(
            transfer_uri=self.transfer_uri)
        self.assertEqual(result, self.subclient)
