import os
import tempfile
import time
import shutil


class BaseSubmissionRunner(object):
    def __init__(self, *args, submission_meta=None, **kwargs):
        self.submission_meta = submission_meta

    def run_submission(self): raise NotImplementedError

    def generate_tmp_dir(self, prefix=None, use_time_suffix=True):
        if prefix is None: prefix = 'sr'
        time_suffix = ''
        if use_time_suffix: time_suffix = '.' + str(time.time())
        scratch_dir = self.submission.get('scratch_dir') or tempfile.mkdtemp()
        tmp_dir = os.path.join(scratch_dir, prefix + time_suffix)
        os.makedirs(tmp_dir)
        return tmp_dir

    def move_to_outputs(self, src=None, outputs_key=None):
        outputs_dest = os.path.join(self.submission['outputs_dir'], outputs_key)
        shutil.move(src, outputs_dest)

def run_job_submission(*args, JobSubmissionRunner=None, submission_meta=None,
                       **kwargs):
    return JobSubmissionRunner(*args, submission_meta=submission_meta,
                               **kwargs).run_submission()
