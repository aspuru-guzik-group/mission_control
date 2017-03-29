import json
import os

from mc.a2g2.conformer_generators.rdkit_conformer_generator\
        .rdkit_conformer_generator import (
            RDKitConformerGenerator as ConformerGenerator)


class SubmissionRunner(object):
    def __init__(self, submission=None):
        self.submission = submission

    def run_submission(self):
        self.generate_conformers(
            output_dir=self.make_output_dir(),
            confgen_params=self.get_confgen_params(),
        )

    def make_output_dir(self):
        output_dir = os.path.join(self.submission['outputs_dir'],
                                  'confgen_outputs')
        os.makedirs(output_dir)
        return output_dir

    def get_confgen_params(self):
        job = self.load_job()
        return job['job_spec']['job_params']['confgen_params']

    def load_job(self):
        job_json_path = os.path.join(self.submission['dir'], 'job.json')
        return json.load(open(job_json_path))

    def generate_conformers(self, confgen_params=None):
        generator = ConformerGenerator(**confgen_params)
        generator.generate_conformers()
