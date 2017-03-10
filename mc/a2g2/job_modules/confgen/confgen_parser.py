import glob
import json
import os


class ConfgenParser(object):
    def parse_completed_confgen_dir(self, completed_confgen_dir=None,
                                    output_dir=None):
        calc_chemthing = self.extract_calc_chemthing(
            completed_confgen_dir=completed_confgen_dir)
        conformer_chemthings = self.extract_conformer_chemthings(
            completed_confgen_dir=completed_confgen_dir,
            calc_chemthing=calc_chemthing)
        self.write_chemthings_bulk_file(
            chemthings=[calc_chemthing, *conformer_chemthings],
            target_path=os.path.join(output_dir, 'chemthings.bulk')
        )
        return output_dir

    def extract_calc_chemthing(self, completed_confgen_dir=None):
        job = self.parse_job_file(completed_confgen_dir=completed_confgen_dir)
        calc_chemthing = {
            'uri': 'a2g2:calculation:some_uri',
            'types': ['a2g2:calculation', 'a2g2:calculation:confgen'],
            'props': {
                'a2g2:confgen:parameters': job['job_spec']['parameters']
            },
            'precursors': [job['job_spec']['precursors']]
        }
        return calc_chemthing

    def parse_job_file(self, completed_confgen_dir=None):
        job_file_path = os.path.join(completed_confgen_dir, 'job.json')
        with open(job_file_path, 'r') as job_file: job = json.load(job_file)
        return job

    def extract_conformer_chemthings(self, completed_confgen_dir=None,
                                     calc_chemthing=None):
        conformer_chemthings = []
        xyz_paths = glob.glob(completed_confgen_dir + '/conformers/*')
        for i, xyz_path in enumerate(xyz_paths):
            with open(xyz_path, 'r') as xyz_file: xyz = xyz_file.read()
            conformer_chemthing = {
                'uri': 'a2g2:molecule3d:uri_%s' % i,
                'types': ['a2g2:molecule', 'a2g2:molecule3d'],
                'props': {
                    'a2g2:xyz': xyz
                },
                'precursors': calc_chemthing['uri'],
            }
            conformer_chemthings.append(conformer_chemthing)
        return conformer_chemthings

    def write_chemthings_bulk_file(self, chemthings=None, target_path=None):
        with open(target_path, 'w') as bulk_file:
            for chemthing in chemthings:
                bulk_file.write(json.dumps(chemthing))
