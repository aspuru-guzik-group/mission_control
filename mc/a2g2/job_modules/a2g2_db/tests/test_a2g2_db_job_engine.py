import unittest
from unittest.mock import call, MagicMock, patch

from ...utils import test_utils as job_module_test_utils
from .. import a2g2_db_job_engine


class BaseTestCase(unittest.TestCase,
                   job_module_test_utils.JobEngineTestCaseMixin):
    def setUp(self):
        unittest.TestCase.setUp(self)
        job_module_test_utils.JobEngineTestCaseMixin.setUp(self)
        self.cfg = MagicMock()

    def generate_engine(self):
        return a2g2_db_job_engine.generate_job_engine()

class UploadBulkFilesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job['job_spec']['command'] = 'upload_bulk_files'

    def decorate_patchers(self):
        self.patchers['a2g2_db_uploader'] = patch.object(
            a2g2_db_job_engine.a2g2_db_uploader, 'generate_uploader')

    def test_calls_upload_chemthings(self):
        self.execute_job(job=self.job, cfg=self.cfg)
        expected_uploader = a2g2_db_job_engine.a2g2_db_uploader.\
                generate_uploader.return_value
        self.assertEqual(
            expected_uploader.upload_bulk_files.call_args,
            call(bulk_file_dir=self.job['job_spec']['kwargs']['bulk_file_dir'],
                 cfg=self.cfg)
        )
