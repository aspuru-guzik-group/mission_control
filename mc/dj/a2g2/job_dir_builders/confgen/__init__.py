import json
import tempfile
import textwrap


class ConfgenJobDirBuilder(object):
    @classmethod
    def build_odyssey_dir(cls, job=None, odyssey_dir_builder=None,
                          output_dir=None):
        if not output_dir: output_dir = tempfile.mkdtemp(prefix='confgen.ody')
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
                source /n/aagfs01/software/envs/a2g2_env/bin/activate
                python -m a2g2.conformer_generator \\
                  --input=./conformer_generator.in.json
                '''),
            'templates': {
                'specs': [
                    {
                        'target': 'confgen.in.json',
                        'content': json.dumps({
                            'smiles': job['spec']['smiles']
                        }, indent=2)
                    }
                ],
            }
        }
        return dir_spec

