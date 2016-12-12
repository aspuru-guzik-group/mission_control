import difflib
import filecmp
import json
import os
import unittest

from ...odyssey.odyssey_dir_builder import OdysseyDirBuilder

TestCase = unittest.TestCase


class BaseTestCaseMixin(object):
    def setUp(self):
        self.builder = OdysseyDirBuilder()

class SnapshotDirTestCaseMixin(BaseTestCaseMixin):
    def setUp(self):
        if not hasattr(self, 'snapshot_name'):
            raise Exception("Define 'snapshot_name' in TestCase class.")
        super().setUp()
        self.snapshot_dir = os.path.join(os.path.dirname(__file__),
                                         self.snapshot_name)
        self.job_spec = self.load_job_spec()
        self.expected_output_dir = self.get_expected_output_dir()

    def load_job_spec(self):
        return json.load(open(os.path.join(self.snapshot_dir, 'job_spec.json')))

    def get_expected_output_dir(self):
        return os.path.join(self.snapshot_dir, 'expected_dir')

    def test_builds_expected_job(self):
        output_dir = self.builder.build_job(job_spec=self.job_spec)
        self.compare_dirs(output_dir, self.expected_output_dir)

    def compare_dirs(self, left, right):
        cmp_result = filecmp.dircmp(left, right)
        self.assertEqual(set(cmp_result.left_only), set())
        self.assertEqual(set(cmp_result.right_only), set())
        for differing_file in cmp_result.diff_files:
            fail_msg = "Files differ:\n%s" % self.get_diff(
                left_file=os.path.join(left, differing_file),
                right_file=os.path.join(right, differing_file))
            self.fail(fail_msg)

    def get_diff(self, left_file=None, right_file=None):
        diff_lines = difflib.unified_diff(
            open(left_file).readlines(), open(right_file).readlines(),
            fromfile=left_file, tofile=right_file)
        return "".join([line for line in diff_lines])

class rdkit_conformer_generator_TestCase(SnapshotDirTestCaseMixin, TestCase):
    snapshot_name = 'rdkit_conformer_generator'

if __name__ == '__main__':
    unittest.main()
