import json
import textwrap
import unittest
from unittest.mock import MagicMock

from .a2g2_job_engine_dir_builder import A2G2JobEngineDirBuilder


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = {'some_key': 'some_value'}
        self.cfg = {'A2G2_JOB_ENGINE_MODULE': 'some_module'}

class TestBuildOdysseyDir(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_passes_expected_dir_spec(self):
        mock_odyssey_builder = MagicMock()
        mock_output_dir = MagicMock()
        A2G2JobEngineDirBuilder.build_odyssey_dir(
            job=self.job, odyssey_dir_builder=mock_odyssey_builder,
            output_dir=mock_output_dir, cfg=self.cfg)
        expected_dir_spec = {
            'job_script_body': textwrap.dedent(
                '''
                set -o errexit
                echo "starting, $(date)"
                source {conda_env_root}/bin/activate {conda_env_root}
                SCRATCH_DIR="/scratch/a2g2_job_engine_job.$$"
                mkdir -p $SCRATCH_DIR
                python -m {a2g2_job_engine_module} \\
                    --job="./{job_file_name}" \\
                    --cfg="./{cfg_file_name}" \\
                    --output_dir="$SCRATCH_DIR" \\
                    --ctx_dir="$(pwd)"
                cp -r $SCRATCH_DIR ./output
                echo "finished, $(date)"
                ''').format(
                    conda_env_root='/n/aagfs01/software/conda_envs/a2g2_env',
                    a2g2_job_engine_module=self.cfg['A2G2_JOB_ENGINE_MODULE'],
                    job_file_name=A2G2JobEngineDirBuilder.job_file_name,
                    cfg_file_name=A2G2JobEngineDirBuilder.cfg_file_name
                ),
            'templates': {
                'specs': [
                    {
                        'target': A2G2JobEngineDirBuilder.job_file_name,
                        'content': json.dumps(self.job, indent=2),
                        'no_render': True,
                    },
                    {
                        'target': A2G2JobEngineDirBuilder.cfg_file_name,
                        'content': json.dumps(self.cfg, indent=2),
                        'no_render': True,
                    }
                ],
            }
        }
        call_kwargs = mock_odyssey_builder.build_dir.call_args[1]
        self.assertEqual(call_kwargs['output_dir'], mock_output_dir)
        for k, v in expected_dir_spec.items():
            self.assertEqual(call_kwargs['dir_spec'][k], v) 
