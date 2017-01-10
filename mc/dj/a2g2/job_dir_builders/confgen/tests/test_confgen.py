import json
import os
import tempfile
import textwrap
import unittest
from unittest.mock import MagicMock

from .. import ConfgenJobDirBuilder
from ...odyssey import OdysseyJobDirBuilder
from ... import test_utils as _test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = {'spec': {'smiles': 'some smiles'}}

class TestBuildOdysseyDir(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_passes_expected_dir_spec(self):
        mock_odyssey_builder = MagicMock()
        mock_output_dir = MagicMock()
        ConfgenJobDirBuilder.build_odyssey_dir(
            job=self.job, odyssey_dir_builder=mock_odyssey_builder,
            output_dir=mock_output_dir)
        expected_dir_spec = {
            'modules': ['conda'],
            'job_script_body': textwrap.dedent(
                '''
                source /n/aagfs01/software/envs/a2g2_env/bin/activate
                python -m a2g2.conformer_generator \\
                  --input=./conformer_generator.in.json
                '''),
            'templates': {
                'specs': [
                    {'target': 'confgen.in.json',
                     'content': json.dumps({
                         'smiles': self.job['spec']['smiles']
                     }, indent=2)}
                ]
            }
        }
        call_kwargs = mock_odyssey_builder.build_dir.call_args[1]
        self.assertEqual(call_kwargs['output_dir'], mock_output_dir)
        for k, v in expected_dir_spec.items():
            self.assertEqual(call_kwargs['dir_spec'][k], v) 

    def test_generates_expected_dir(self):
        output_dir = tempfile.mkdtemp(prefix='confgen.ody.')
        ConfgenJobDirBuilder.build_odyssey_dir(
            job=self.job, odyssey_dir_builder=OdysseyJobDirBuilder,
            output_dir=output_dir)
        this_dir = os.path.dirname(__file__)
        snapshot_dir = os.path.join(this_dir, 'snapshots', 'odyssey',
                                    'expected_dir')
        _test_utils.assert_dirs_equal(self, output_dir, snapshot_dir)
