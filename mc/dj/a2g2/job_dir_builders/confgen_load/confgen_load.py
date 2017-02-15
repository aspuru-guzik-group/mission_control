import json
import shlex
import tempfile
import textwrap

from ..odyssey import OdysseyJobDirBuilder


class ConfgenLoadJobDirBuilder(object):
    @classmethod
    def build_odyssey_dir(cls, job=None, odyssey_dir_builder=None,
                          output_dir=None, a2g2_client_cfg_json=None):
        if not output_dir: output_dir = tempfile.mkdtemp(
            prefix='confgen_load.ody')
        if not odyssey_dir_builder: odyssey_dir_builder = OdysseyJobDirBuilder
        dir_spec = cls.generate_odyssey_dir_spec(
            job=job,
            a2g2_client_cfg_json=a2g2_client_cfg_json
        )
        dir_path = odyssey_dir_builder.build_dir(
            dir_spec=dir_spec,
            output_dir=output_dir
        )
        return dir_path

    @classmethod
    def generate_odyssey_dir_spec(cls, job=None, a2g2_client_cfg_json=None):
        dir_spec = {
            'modules': ['conda'],
            'job_script_body': textwrap.dedent(
                '''
                echo "starting, $(date)"
                source activate /n/aagfs01/software/conda_envs/a2g2_env
                export A2G2_CLIENT_CFG_JSON="{a2g2_client_cfg_json}"
                python -m a2g2.job_engines.confgen_load_job_engine \\
                    --job="./job.json" \\
                    --cfg="./cfg.json" \\
                echo "finished, $(date)"
                '''
            ).format(a2g2_client_cfg_json=shlex.quote(a2g2_client_cfg_json)),
            'templates': {
                'specs': [
                    {
                        'target': 'job.json',
                        'content': json.dumps(job, indent=2),
                        'no_render': True,
                    }
                ]
            }
        }
        return dir_spec

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Build confgen load dir.")
    parser.add_argument('--job', type=lambda s: json.load(open(s)))
    parser.add_argument('--output_dir', type=str)
    args = vars(parser.parse_args())
    ConfgenLoadJobDirBuilder.build_odyssey_dir(**args)
    
