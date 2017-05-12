import hashlib
import json
import os

from .. import constants as qchem_constants


class QChemParser(object):
    def __init__(self, *args, input_dir=None, output_dir=None,
                 parsing_params=None, **kwargs):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.parsing_params = parsing_params or {}
        self.completed_qchem_workdir = self.get_completed_qchem_workdir()
        self.qchem_output_file_path = os.path.join(
            self.completed_qchem_workdir,
            qchem_constants.QCHEM_OUTPUT_FILE_NAME)

    def get_completed_qchem_workdir(self):
        submission = self.parse_submission()
        return os.path.join(self.input_dir, submission['outputs_dir'],
                            qchem_constants.QCHEM_OUTPUTS_KEY)

    def parse_submission(self):
        submission_file_path = os.path.join(self.input_dir, 'submission.json')
        return json.load(open(submission_file_path))

    def parse(self):
        calc_chemthing = self.extract_calc_chemthing()
        self.write_chemthings_bulk_file(
            chemthings=[calc_chemthing],
            target_path=os.path.join(self.output_dir, 'qchem.chemthings.bulk')
        )
        return self.output_dir

    def extract_calc_chemthing(self):
        calc_chemthing = {
            'types': {'a2g2:type:calc': True},
            'keys': {self.generate_calculation_key(): True},
            'props': {
                'a2g2:calculation:artifact': self.parsing_params['artifact'],
                'a2g2:calculation:qchem:rem_items': self.extract_rem_items()
            },
            'precursors': self.parsing_params.get('precursors')
        }
        return calc_chemthing

    def generate_calculation_key(self):
        return 'a2g2:calculation:{hash}'.format(
            hash=self.get_hash_for_file(file_path=self.qchem_output_file_path))

    def get_hash_for_file(self, file_path=None):
        BUF_SIZE = 65536
        hasher = hashlib.sha256()
        with open(file_path, 'rb') as f:
            while True:
                data = f.read(BUF_SIZE)
                if not data: break
                hasher.update(data)
        return hasher.hexdigest()

    def extract_rem_items(self):
        rem_items = []
        with open(self.qchem_output_file_path) as f:
            in_rem_section = False
            for line in f:
                line = line.strip()
                if not line or line.startswith('!'): continue
                if line.startswith('$rem'):
                    in_rem_section = True
                    continue
                if in_rem_section:
                    if line.startswith('$end'): break
                    rem_item_parts = line.split()
                    rem_item = {
                        'variable': rem_item_parts[0],
                        'value': rem_item_parts[1],
                    }
                    if len(rem_item_parts) == 3:
                        rem_item['comment'] = rem_item_parts[2]
                    rem_items.append(rem_item)
        return rem_items

    def write_chemthings_bulk_file(self, chemthings=None, target_path=None):
        with open(target_path, 'w') as bulk_file:
            for chemthing in chemthings:
                bulk_file.write(json.dumps(chemthing) + "\n")

def parse_qchem_dir(*args, transform_params=None, **kwargs):
    parser = QChemParser(*args, parsing_params=transform_params, **kwargs)
    return parser.parse()
