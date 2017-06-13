import os
import tempfile
import shutil


class BaseJobSubmissionRunner(object):
    def __init__(self, *args, submission_meta=None, **kwargs):
        self.submission_meta = submission_meta

    @classmethod
    def run_job_submission(cls, *args, submission_meta=None, **kwargs):
        runner = cls(*args, submission_meta=submission_meta, **kwargs)
        return runner._run_job_submission()

    def _run_job_submission(self): raise NotImplementedError

    def generate_tmp_dir(self, prefix=None, suffix=None):
        scratch_dir = self.submission_meta.get('scratch_dir')
        if not scratch_dir:
            scratch_dir = tempfile.mkdtemp(prefix='MC_SCRATCH.',
                                           dir=self.submission_meta['dir'])
        return tempfile.mkdtemp(prefix=prefix, suffix=suffix, dir=scratch_dir)

    def move_to_outputs(self, src=None, outputs_key=None):
        outputs_dir = os.path.join(
            self.submission_meta['dir'],
            self.submission_meta.get('io_dir_names', {}).get(
                'outputs', 'outputs')
        )
        outputs_dest = os.path.join(outputs_dir, outputs_key)
        shutil.move(src, outputs_dest)

    @property
    def job_params(self):
        if not hasattr(self, '_job_params'):
            self.job_params = \
                    self.submission_meta['job']['job_spec']['job_params']
        return self._job_params

    @job_params.setter
    def job_params(self, value): self._job_params = value
