import os
import tempfile
import time
import shutil


class BaseSubmissionRunner(object):
    def __init__(self, job=None, cfg=None, submission=None):
        self.job = job
        self.cfg = cfg
        self.submission = submission

    def run_submission(self): raise NotImplementedError

    def generate_tmp_dir(self, prefix=None, use_time_suffix=True):
        if prefix is None: prefix = ''
        time_suffix = ''
        if use_time_suffix: time_suffix = str(time.time())
        scratch_dir = self.submission.get('scratch_dir')
        if scratch_dir:
            tmp_dir = os.path.join(scratch_dir, prefix)
        else:
            tmp_dir = tempfile.mkdtemp()
        tmp_dir = tmp_dir + time_suffix
        os.makedirs(tmp_dir)
        return tmp_dir

    def move_to_outputs(self, src=None, outputs_key=None):
        outputs_dest = os.path.join(self.submission['outputs_dir'], outputs_key)
        shutil.move(src, outputs_dest)

def run_job_submission(*args, JobSubmissionRunner=None, job=None, cfg=None,
                       submission=None, **kwargs):
    return JobSubmissionRunner(*args, job=job, cfg=cfg, submission=submission,
                               **kwargs).run_submission()
