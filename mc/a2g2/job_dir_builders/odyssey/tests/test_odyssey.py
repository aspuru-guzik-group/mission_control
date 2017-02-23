import json
import os
import tempfile
import unittest

from .. import OdysseyJobDirBuilder
from ... import test_utils as _test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self): pass

class BuildDirTestCase(BaseTestCase):
    def test_builds_expected_dir(self):
        this_dir = os.path.dirname(__file__)
        snapshot_dir = os.path.join(this_dir, 'snapshots', 'spec_1')
        dir_spec = json.load(open(os.path.join(snapshot_dir, 'dir_spec.json')))
        output_dir = tempfile.mkdtemp(prefix='ody.test.')
        dir_meta = OdysseyJobDirBuilder.build_dir(dir_spec=dir_spec,
                                                output_dir=output_dir)
        expected_dir_meta = {
            'dir': output_dir, 
            'checkpoint_files': OdysseyJobDirBuilder.checkpoint_files,
            'output_files': OdysseyJobDirBuilder.output_files,
            'entrypoint': 'job.sh',
        }
        self.assertEqual(dir_meta, expected_dir_meta)
        expected_dir = os.path.join(snapshot_dir, 'expected_output')
        _test_utils.assert_dirs_equal(self, output_dir, expected_dir)

if __name__ == '__main__':
    unittest.main()
