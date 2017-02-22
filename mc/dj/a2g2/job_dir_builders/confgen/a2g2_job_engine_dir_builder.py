import json
import tempfile
import textwrap

from ..odyssey import OdysseyJobDirBuilder


class A2G2JobEngineDirBuilder(object):
    job_file_name = 'job.json'
    cfg_file_name = 'engine_cfg.json'

    @classmethod
    def build_odyssey_dir(cls, job=None, odyssey_dir_builder=None,
                          output_dir=None, cfg=None):
        if not output_dir: output_dir = tempfile.mkdtemp(prefix='a2g2jedb.ody')
        if not odyssey_dir_builder: odyssey_dir_builder = OdysseyJobDirBuilder
        dir_spec = cls.generate_odyssey_dir_spec(job=job, cfg=cfg)
        dir_meta = odyssey_dir_builder.build_dir(dir_spec=dir_spec,
                                                 output_dir=output_dir)
        return dir_meta

    @classmethod
    def generate_odyssey_dir_spec(cls, job=None, cfg=None):
        dir_spec = {
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
                    --output_dir="$SCRATCH_DIR"
                cp -r $SCRATCH_DIR ./output
                echo "finished, $(date)"
                ''').format(
                    conda_env_root='/n/aagfs01/software/conda_envs/a2g2_env',
                    a2g2_job_engine_module='a2g2.job_engines.a2g2_job_engine',
                    job_file_name=cls.job_file_name,
                    cfg_file_name=cls.cfg_file_name
                ),
            'templates': {
                'specs': [
                    {
                        'target': cls.job_file_name,
                        'content': json.dumps(job, indent=2)
                    },
                    {
                        'target': cls.cfg_file_name,
                        'content': json.dumps(cfg, indent=2)
                    }
                ],
            }
        }
        return dir_spec

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Build a2g2_job_engine dir.")
    parser.add_argument('--job', type=lambda s: json.load(open(s)))
    parser.add_argument('--output_dir', type=str)
    args = vars(parser.parse_args())
    A2G2JobEngineDirBuilder.build_odyssey_dir(**args)
    
