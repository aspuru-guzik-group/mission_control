import json
import tempfile
import textwrap

from ..odyssey import OdysseyJobDirBuilder


class ConfgenJobDirBuilder(object):
    @classmethod
    def build_odyssey_dir(cls, job=None, odyssey_dir_builder=None,
                          output_dir=None):
        if not output_dir: output_dir = tempfile.mkdtemp(prefix='confgen.ody')
        if not odyssey_dir_builder: odyssey_dir_builder = OdysseyJobDirBuilder
        dir_spec = cls.generate_odyssey_dir_spec(job=job)
        dir_path = odyssey_dir_builder.build_dir(dir_spec=dir_spec,
                                                 output_dir=output_dir)
        return dir_path

    @classmethod
    def generate_odyssey_dir_spec(cls, job=None):
        dir_spec = {
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
                    {
                        'target': 'confgen.params.json',
                        'content': json.dumps({
                            'smiles': job['job_spec']['confgen']['smiles']
                        }, indent=2)
                    }
                ],
            }
        }
        return dir_spec

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Build confgen dir.")
    parser.add_argument('--job', type=lambda s: json.load(open(s)))
    parser.add_argument('--output_dir', type=str)
    args = vars(parser.parse_args())
    ConfgenJobDirBuilder.build_odyssey_dir(**args)
    
