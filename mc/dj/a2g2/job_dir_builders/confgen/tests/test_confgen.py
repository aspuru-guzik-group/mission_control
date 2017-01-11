import json
import textwrap
import unittest
from unittest.mock import MagicMock

from ..confgen import ConfgenJobDirBuilder


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
                echo "starting, $(date)"
                source activate /n/aagfs01/software/conda_envs/a2g2_env
                SCRATCH_DIR="/scratch/conformers.$$"
                mkdir -p $SCRATCH_DIR
                python -m a2g2_utils.conformer_generators.rdkit_conformer_generator.cmd \\
                    --params_file="./confgen.params.json" \\
                    --output_dir="$SCRATCH_DIR"
                cp -r $SCRATCH_DIR ./conformers
                echo "finished, $(date)"
                '''),
            'templates': {
                'specs': [
                    {'target': 'confgen.params.json',
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
