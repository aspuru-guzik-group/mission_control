import os
import tempfile
import time
import shutil


class BaseJobSubmissionRunner(object):
    def __init__(self, *args, submission_meta=None, **kwargs):
        self.submission_meta = submission_meta

    @classmethod
    def run_job_submission(cls, *args, submission_meta=None, **kwargs):
        runner = cls(*args, submission_meta=submission_meta, **kwargs)
        return runner._run_job_submission()

    def _run_job_submission(self): raise NotImplementedError

    def generate_tmp_dir(self, prefix=None, use_time_suffix=True):
        if prefix is None: prefix = 'sr'
        time_suffix = ''
        if use_time_suffix: time_suffix = '.' + str(time.time())
        scratch_dir = self.submission_meta.get('scratch_dir') \
                or tempfile.mkdtemp()
        tmp_dir = os.path.join(scratch_dir, prefix + time_suffix)
        os.makedirs(tmp_dir)
        return tmp_dir

    def move_to_outputs(self, src=None, outputs_key=None):
        outputs_dir = os.path.join(
            self.submission_meta['dir'],
            self.submission_meta.get('io_dir_names', {}).get(
                'outputs', 'outputs')
        )
        outputs_dest = os.path.join(outputs_dir, outputs_key)
        shutil.move(src, outputs_dest)
