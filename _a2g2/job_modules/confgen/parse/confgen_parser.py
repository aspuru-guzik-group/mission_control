import glob
import json
import os
from uuid import uuid4

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
                            confgen_constants.CONFGEN_OUTPUTS_KEY)

    def parse_submission(self):
        submission_file_path = os.path.join(self.input_dir, 'submission.json')
        return json.load(open(submission_file_path))

    def parse_confgen_dir(self):
        calc_chemthing = self.extract_calc_chemthing()
        conformer_chemthings = self.extract_conformer_chemthings(
            calc_chemthing=calc_chemthing)
        self.write_chemthings_bulk_file(
            chemthings=[calc_chemthing, *conformer_chemthings],
            target_path=os.path.join(self.output_dir, 'confgen.chemthings.bulk')
        )
        return self.output_dir

    def extract_calc_chemthing(self):
        calc_chemthing = {
            'uuid': str(uuid4()),
            'types': {
                'a2g2:type:calc:confgen': True,
            },
            'props': {
                'a2g2:prop:confgen:parameters': self.parse_confgen_params(),
            },
            'precursors': self.parsing_params.get('precursors', {})
        }
        return calc_chemthing

    def parse_confgen_params(self):
        file_path = os.path.join(self.completed_confgen_workdir, 'outputs',
                                 confgen_constants.CONFGEN_PARAMS_FILENAME)
        return json.load(open(file_path))

    def extract_conformer_chemthings(self, calc_chemthing=None):
        conformer_chemthings = []
        xyz_dir = os.path.join(self.completed_confgen_workdir, 'outputs',
                               'conformers')
        xyz_paths = glob.glob(xyz_dir + '/*.xyz')
        for i, xyz_path in enumerate(xyz_paths):
            parsed_xyz = self.parse_xyz(xyz=open(xyz_path).read())
            conformer_chemthing = {
                'uuid': str(uuid4()),
                'types': {'a2g2:type:mol3d': True},
                'props': {
                    'a2g2:prop:atoms': parsed_xyz['atoms'],
                    'a2g2:prop:confgen:comment': parsed_xyz['comment'],
                },
                'precursors': {calc_chemthing['uuid']: True},
            }
            conformer_chemthings.append(conformer_chemthing)
        return conformer_chemthings

    def parse_xyz(self, xyz=None):
        xyz_lines = xyz.split("\n")
        parsed_xyz = {
            'num_atoms': int(xyz_lines[0]),
            'comment': xyz_lines[1],
            'atoms': [self.parse_atom_line(atom_line)
                      for atom_line in xyz_lines[2:]]
        }
        return parsed_xyz

    def parse_atom_line(self, atom_line=None):
        line_parts = atom_line.split()
        element = line_parts[0]
        x, y, z = [float("%.4f" % float(line_part))
                   for line_part in line_parts[1:]]
        return {'element': element, 'x': x, 'y': y, 'z': z}

    def write_chemthings_bulk_file(self, chemthings=None, target_path=None):
        with open(target_path, 'w') as bulk_file:
            for chemthing in chemthings:
                bulk_file.write(json.dumps(chemthing) + "\n")

def parse_confgen_dir(*args, transform_params=None, **kwargs):
    parser = ConfgenParser(*args, parsing_params=transform_params, **kwargs)
    return parser.parse_confgen_dir()
