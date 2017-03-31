import os
import shutil

from .. import base_submission_runner


class BaseParseSubmissionRunner(base_submission_runner.BaseSubmissionRunner):
    def __init__(self, *args, parse_dir_fn=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.parse_dir_fn = parse_dir_fn or self.parse_dir_fn

    def parse_dir_fn(self, input_dir=None, output_dir=None):
        raise NotImplementedError

    def run_submission(self):
        outdir = self.generate_scratchdir()
        self.parse_dir_fn(input_dir=self.get_input_dir(), output_dir=outdir)
        outputs_dest = os.path.join(self.submission['outputs_dir'], 'parse')
        shutil.move(outdir, outputs_dest)

    def get_input_dir(self):
        return os.path.join(self.submission['inputs_dir'], 'dir_to_parse')

def run_job_submission(*args, parse_dir_fn=None, **kwargs):
    return base_submission_runner.run_job_submission(
        *args,
        JobSubmissionRunner=BaseParseSubmissionRunner,
        parse_dir_fn=parse_dir_fn,
        **kwargs)
                    
