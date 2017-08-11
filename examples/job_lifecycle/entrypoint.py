import tempfile

from a2g2_v2.space_shuttle.space_shuttle import SpaceShuttle


class Entrypoint(object):
    def run(self):
        self.shuttle = SpaceShuttle(ensure_db=False, ensure_job_dirs=False)
        job_dict = {
            'key': 'some-job-key',
            'module': 'a2g2_v2.confgen.job_module',
            'params': {
                'smiles': 'C1=CC(=C(C=C1CCN)O)O',
            }
        }
        output_dir = tempfile.mkdtemp(prefix='confgen.job_dir.')
        self.shuttle.run_command(
            'build_job_dir',
            job_dict=job_dict,
            output_dir=output_dir
        )
        return output_dir


if __name__ == '__main__':
    Entrypoint().run()
