import json
import textwrap
import unittest
from unittest.mock import MagicMock

from .confgen_load import ConfgenLoadJobDirBuilder


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = {
            'job_spec': {
                'confgen_load': {
                }
            }
        }

class TestBuildOdysseyDir(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_passes_expected_dir_spec(self):
        mock_odyssey_builder = MagicMock()
        mock_output_dir = MagicMock()
        a2g2_client_cfg_json = 'mock_cfg'
        ConfgenLoadJobDirBuilder.build_odyssey_dir(
            job=self.job, odyssey_dir_builder=mock_odyssey_builder,
            output_dir=mock_output_dir,
            a2g2_client_cfg_json=a2g2_client_cfg_json
        )
        expected_dir_spec = {
            'modules': ['conda'],
            'job_script_body': textwrap.dedent(
                '''
                echo "starting, $(date)"
                source activate /n/aagfs01/software/conda_envs/a2g2_env
                export A2G2_CLIENT_CFG_JSON={a2g2_client_cfg_json}
                python -m a2g2.job_engines.confgen_load_job_engine \\
                    --job="./job.json" \\
                    --cfg="./cfg.json" \\
                echo "finished, $(date)"
                '''
            ).format(a2g2_client_cfg_json=a2g2_client_cfg_json),
            'templates': {
                'specs': [
                    {
                        'target': 'job.json',
                        'content': json.dumps(self.job, indent=2)
                    }
                ]
            }
        }
        call_kwargs = mock_odyssey_builder.build_dir.call_args[1]
        self.assertEqual(call_kwargs['output_dir'], mock_output_dir)
        for k, v in expected_dir_spec.items():
            self.assertEqual(call_kwargs['dir_spec'][k], v) 
