import json
import os
import re
from uuid import uuid4

from .. import constants as qchem_constants


class QChemComputationParser(object):
    def __init__(self, *args, parsing_params=None, **kwargs):
        self.parsing_params = parsing_params

    def extract_computation_meta(self, job_dir=None):
        computation_meta = {
            'uuid': str(uuid4()),
            'command_meta': self.extract_command_meta(
                qchem_output_file=self.get_qchem_output_file(job_dir=job_dir)),
            'execution_meta': self.extract_execution_meta(job_dir=job_dir)
        }
        return computation_meta

    def get_qchem_output_file(self, job_dir=None):
        completed_qchem_work_dir = self.get_completed_qchem_work_dir(
            job_dir=job_dir)
        self.qchem_output_file_path = os.path.join(
            completed_qchem_work_dir,
            qchem_constants.QCHEM_OUTPUT_FILE_NAME)
        return open(self.qchem_output_file_path)

    def get_completed_qchem_work_dir(self, job_dir=None):
        submission = self.parse_submission(job_dir=job_dir)
        return os.path.join(job_dir, submission['outputs_dir'],
                            qchem_constants.QCHEM_OUTPUTS_KEY)

    def parse_submission(self, job_dir=None):
        submission_file_path = os.path.join(job_dir, 'submission.json')
        return json.load(open(submission_file_path))

    def extract_command_meta(self, qchem_output_file=None):
        command_meta = {
            'executable': self.extract_executable_meta(
                qchem_output_file=qchem_output_file),
            'qchem_params': {
                rem_item['variable']: rem_item['value']
                for rem_item in self.extract_rem_items(
                    qchem_output_file=qchem_output_file)
            }
        }
        return command_meta

    def extract_executable_meta(self, qchem_output_file=None):
        executable_meta = {'name': 'Q-Chem'}
        for line in qchem_output_file:
            line = line.strip()
            m = re.match(r'Q-Chem (\d+\.\d+\.\d+)', line)
            if m:
                executable_meta['version'] = m.group(1)
                break
        return executable_meta

    def extract_rem_items(self, qchem_output_file=None):
        rem_items = []
        in_rem_section = False
        for line in qchem_output_file:
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

    def extract_execution_meta(self, job_dir=None):
        submission = self.parse_submission(job_dir=job_dir)
        if 'execution_meta_name' not in submission: return {}
        execution_meta_path = os.path.join(job_dir, submission['outputs_dir'],
                                           submission['execution_meta_name'])
        return json.load(open(execution_meta_path))
