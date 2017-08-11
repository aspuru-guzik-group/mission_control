import unittest

from mc.examples import testing_utils
from . import entrypoint


class TestEntrypoint(unittest.TestCase):
    def test_entrypoint(self):
        output_dir = entrypoint.Entrypoint().run()
        output_dir_items = testing_utils.get_dir_items(output_dir)
        expected_dir_items = [
            'A2G2__JOB_KEY', 'A2G2__JOB_META.json', 'JOBMAN__JOB_SPEC.json',
            'entrypoint.sh', 'work_dir', 'work_dir/confgen.in.json',
            'work_dir/entrypoint.sh'
        ]
        self.assertEqual(output_dir_items, expected_dir_items)
