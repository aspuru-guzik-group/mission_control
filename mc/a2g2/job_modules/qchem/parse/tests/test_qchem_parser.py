import json
import gzip
import hashlib
import os
import tempfile
import unittest

from mc.mc_utils import test_utils as mc_test_utils

from ... import constants as qchem_constants
from .. import qchem_parser


class ParseQChemDirTestCase(unittest.TestCase):
    def setUp(self):
        self.parsing_params = {
            'artifact': {'artifact1': 'artifact1_meta'},
            'precursors': {'precursor1': True}
        }
        self.output_dir = tempfile.mkdtemp()
        self.completed_qchem_dir_outputs_name = 'outputs'
        self.completed_qchem_dir  = self.generate_completed_qchem_dir()
        self.expected_output_dir = self.generate_expected_output_dir(
            completed_qchem_dir=self.completed_qchem_dir)

    def generate_completed_qchem_dir(self):
        completed_qchem_dir = tempfile.mkdtemp()
        outputs_dir = os.path.join(completed_qchem_dir,
                                   self.completed_qchem_dir_outputs_name)
        os.makedirs(outputs_dir)
        self.generate_submission_file(parent_dir=completed_qchem_dir,
                                      outputs_dir=outputs_dir)
        self.generate_qchem_output_dir(outputs_dir=outputs_dir)
        return completed_qchem_dir

    def generate_submission_file(self, parent_dir=None, outputs_dir=None):
        submission = {'outputs_dir': outputs_dir}
        with open(os.path.join(parent_dir, 'submission.json'), 'w') as f:
            json.dump(submission, f)

    def generate_qchem_output_dir(self, outputs_dir=None):
        qchem_output_dir = os.path.join(outputs_dir,
                                        qchem_constants.QCHEM_OUTPUTS_KEY)
        os.makedirs(qchem_output_dir)
        qchem_output_file_path = os.path.join(
            qchem_output_dir, qchem_constants.QCHEM_OUTPUT_FILE_NAME)
        example_qchem_output_path = os.path.join(
            os.path.dirname(__file__), 'example_qchem.out.gz')
        open(qchem_output_file_path, 'wb').write(
            gzip.open(example_qchem_output_path, 'rb').read())
        return qchem_output_dir

    def generate_expected_output_dir(self, completed_qchem_dir=None):
        output_dir = tempfile.mkdtemp()
        calc_chemthing = self.generate_expected_calc_chemthing(
            completed_qchem_dir=completed_qchem_dir)
        self.write_chemthings_bulk_file(
            chemthings=[calc_chemthing],
            target_path=os.path.join(output_dir, 'qchem.chemthings.bulk')
        )
        return output_dir

    def generate_expected_calc_chemthing(self, completed_qchem_dir=None):
        expected_key = self.generate_expected_key(
            completed_qchem_dir=completed_qchem_dir)
        calc_chemthing = {
            'types': {'a2g2:type:calc': True},
            'keys': {expected_key: True},
            'props': {
                'a2g2:calculation:artifact': self.parsing_params['artifact'],
                'a2g2:calculation:qchem:rem_items': self.parse_rem_section(
                    completed_qchem_dir=completed_qchem_dir)
            },
            'precursors': self.parsing_params['precursors']
        }
        return calc_chemthing

    def generate_expected_key(self, completed_qchem_dir=None):
        qchem_output_file_path = self.get_qchem_output_file_path(
            completed_qchem_dir=completed_qchem_dir)
        expected_key = "a2g2:calculation:{}".format(self.get_hash_for_file(
            file_path=qchem_output_file_path))
        return expected_key

    def get_qchem_output_file_path(self, completed_qchem_dir=None):
        return os.path.join(completed_qchem_dir,
                            self.completed_qchem_dir_outputs_name,
                            qchem_constants.QCHEM_OUTPUTS_KEY,
                            qchem_constants.QCHEM_OUTPUT_FILE_NAME)

    def get_hash_for_file(self, file_path=None):
        BUF_SIZE = 65536  # lets read stuff in 64kb chunks!
        hasher = hashlib.sha256()
        with open(file_path, 'rb') as f:
            while True:
                data = f.read(BUF_SIZE)
                if not data: break
                hasher.update(data)
        return hasher.hexdigest()

    def parse_rem_section(self, completed_qchem_dir=None):
        qchem_output_file_path = self.get_qchem_output_file_path(
            completed_qchem_dir=completed_qchem_dir)
        rem_items = []
        with open(qchem_output_file_path) as f:
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
                bulk_file.write(json.dumps(chemthing))

    def test_generates_expected_output_dir(self):
        qchem_parser.parse_qchem_dir(input_dir=self.completed_qchem_dir,
                                     output_dir=self.output_dir,
                                     transform_params=self.parsing_params)
        mc_test_utils.assert_dirs_equal(
            test=self,
            left=self.output_dir,
            right=self.expected_output_dir,
            json_patterns=[r'\.json$', r'chemthings\.bulk$']
        )
