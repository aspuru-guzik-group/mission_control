from collections import defaultdict
import json
import os
import tempfile
import unittest
from unittest.mock import call, MagicMock, patch

from .. import a2g2_db_uploader


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.uploader = a2g2_db_uploader.A2G2_DB_Uploader()

class UploadFromBulkFilesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cfg = defaultdict(MagicMock)
        self.cfg['a2g2_client_cfg_json'] = json.dumps({})
        self.chemthings = self.generate_chemthings()
        self.dir_to_upload = self.generate_dir_to_upload(
            chemthings=self.chemthings)
        self.uploader.generate_a2g2_client = MagicMock()
        self.expected_client = self.uploader.generate_a2g2_client.return_value
        self.result = self.uploader.upload_bulk_files(
            bulk_files_dir=self.dir_to_upload, cfg=self.cfg)

    def generate_chemthings(self):
        chemthings = [{'uri': 'chemthing_%s' % i} for i in range(3)]
        return chemthings

    def generate_dir_to_upload(self, chemthings=None):
        _dir = tempfile.mkdtemp()
        chemthings_bulk_path = os.path.join(_dir, 'test.chemthings.bulk')
        with open(chemthings_bulk_path, 'w') as chemthings_bulk_file:
            for chemthing in chemthings:
                chemthings_bulk_file.write(json.dumps(chemthing) + "\n")
        return _dir

    def test_uploads_items_via_client(self):
        expected_call_args_list = [call(chemthing=chemthing)
                                   for chemthing in self.chemthings]
        self.assertEqual(self.expected_client.create_chemthing.call_args_list,
                         expected_call_args_list)

    def test_outputs_upload_results(self):
        expected_result = [self.expected_client.create_chemthing.return_value
                           for chemthing in self.chemthings]
        self.assertEqual(self.result, expected_result)

class GenerateA2G2ClientTestCase(BaseTestCase):
    @patch.object(a2g2_db_uploader, 'a2g2_db_utils')
    def test_dispatches_to_a2g2_db_utils(self, mock_db_utils):
        args = [MagicMock() for i in range(3)]
        kwargs = {str(i): MagicMock() for i in range(3)}
        result = self.uploader.generate_a2g2_client(*args, **kwargs)
        self.assertEqual(mock_db_utils.generate_a2g2_client.call_args,
                         call(*args, **kwargs))
        self.assertEqual(result, mock_db_utils.generate_a2g2_client.return_value)

