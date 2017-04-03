import os
import shutil

from . import base_submission_runner


class BaseTransformSubmissionRunner(base_submission_runner.BaseSubmissionRunner):
    def __init__(self, *args, transform_fn=None, output_key=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.transform_fn = transform_fn or self.transform_fn
        self.output_key = output_key or 'output'

    def run_submission(self):
        scratch_dir = self.generate_scratchdir()
        self.transform_fn(
            input_dir=self.get_input_dir(),
            output_dir=scratch_dir,
            transform_params=self.job['job_spec'].get('job_params'),
            cfg=self.cfg,
            job=self.job,
        )
        outputs_dest = os.path.join(self.submission['outputs_dir'],
                                    self.output_key)
        shutil.move(scratch_dir, outputs_dest)

    def transform_fn(self, input_dir=None, output_dir=None):
        raise NotImplementedError

    def get_input_dir(self):
        return os.path.join(self.submission['inputs_dir'], 'input_dir')

def run_job_submission(*args, transform_fn=None, output_key=None, **kwargs):
    return base_submission_runner.run_job_submission(
        *args,
        JobSubmissionRunner=BaseTransformSubmissionRunner,
        transform_fn=transform_fn,
        output_key=output_key,
        **kwargs)
                    
