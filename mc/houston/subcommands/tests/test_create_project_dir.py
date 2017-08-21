from pathlib import Path
import tempfile
import textwrap
import uuid
import unittest

from mc.houston.tests import utils as _houston_test_utils


class CreateDirTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        houston = _houston_test_utils.generate_test_houston()
        cls.scratch_dir = tempfile.mkdtemp()
        cls.output_dir = str(Path(cls.scratch_dir, str(uuid.uuid4())))
        return houston.run_command(
            'create_project_dir', output_path=cls.output_dir)

    def test_creates_dir(self):
        self.assertTrue(Path(self.output_dir).exists())

    def test_creates_job_dirs_dir(self):
        expected_jobs_dir_path = Path(self.output_dir, 'job_dirs')
        self.assertTrue(expected_jobs_dir_path.exists())

    def test_writes_expected_cfg(self):
        expected_cfg_path = Path(self.output_dir, 'houston_cfg.py')
        self.assertTrue(expected_cfg_path.exists())
        expected_cfg_content = textwrap.dedent(
            '''
            from pathlib import Path

            dir_path = Path(__file__).parent.absolute()

            MC_DB_URI = 'sqlite:///%s' % (dir_path / 'mc.db.sqlite')

            JOB_DIRS_ROOT = str(dir_path / 'job_dirs')
            '''
        ).lstrip()
        with expected_cfg_path.open() as f:
            cfg_content = f.read()
        self.assertEqual(cfg_content, expected_cfg_content)
