from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import default_job_module_loader


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = MagicMock()
        self.cfg = defaultdict(MagicMock)
        self.loader = default_job_module_loader.DefaultJobModuleLoader()

class DefaultJobModuleLoaderCase(BaseTestCase):
    @patch.object(default_job_module_loader, 'importlib')
    def test_imports_expected_module(self, patched_importlib):
        self.loader.get_job_module_name = MagicMock()
        result = self.loader.load_job_module(job=self.job, cfg=self.cfg)
        self.assertEqual(patched_importlib.import_module.call_args,
                         call(self.loader.get_job_module_name.return_value))
        self.assertEqual(result, patched_importlib.import_module.return_value)

    def test_gets_job_module_name_from_job_spec(self):
        result = self.loader.get_job_module_name(job=self.job)
        self.assertEqual(result, self.job['job_spec']['job_type'])
