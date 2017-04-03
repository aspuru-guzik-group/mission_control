import glob
import json
import os

from .. import constants as confgen_constants


class ConfgenParser(object):
    def __init__(self, *args, input_dir=None, output_dir=None,
                 parsing_params=None, **kwargs):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.parsing_params = parsing_params or {}
        self.completed_confgen_workdir = self.get_completed_confgen_workdir()

    def get_completed_confgen_workdir(self):
        submission = self.parse_submission()
        return os.path.join(self.input_dir, submission['outputs_dir'],
                            confgen_constants.OUTPUTS_KEY)

    def parse_submission(self):
        submission_file_path = os.path.join(self.input_dir, 'submission.json')
        return json.load(open(submission_file_path))

    def parse_confgen_dir(self):
        calc_chemthing = self.extract_calc_chemthing()
        conformer_chemthings = self.extract_conformer_chemthings(
            calc_chemthing=calc_chemthing)
        self.write_chemthings_bulk_file(
            chemthings=[calc_chemthing, *conformer_chemthings],
            target_path=os.path.join(self.output_dir, 'chemthings.bulk')
        )
        return self.output_dir

    def extract_calc_chemthing(self):
        confgen_inputs = self.parse_confgen_input_file()
        calc_chemthing = {
            'uri': 'a2g2:calculation:some_uri',
            'types': ['a2g2:calculation', 'a2g2:calculation:confgen'],
            'props': {
                'a2g2:confgen:parameters': confgen_inputs,
            },
            'precursors': self.parsing_params.get('precursors')
        }
        return calc_chemthing

    def parse_confgen_input_file(self):
        file_path = os.path.join(self.completed_confgen_workdir,
                                 confgen_constants.CONFGEN_IN_FILENAME)
        return json.load(open(file_path))

    def extract_conformer_chemthings(self, calc_chemthing=None):
        conformer_chemthings = []
        xyz_paths = glob.glob(self.completed_confgen_workdir + '/conformers/*')
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

def parse_confgen_dir(*args, **kwargs):
    return ConfgenParser(*args, **kwargs).parse_confgen_dir()
