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

    def generate_scratchdir(self):
        scratch_dir = self.submission.get('scratch_dir') or tempfile.mkdtemp()
        workdir = os.path.join(scratch_dir, 'confgen.%s' % time.time())
        os.makedirs(workdir)
        return workdir

    def move_to_outputs(self, src=None, outputs_key=None):
        outputs_dest = os.path.join(self.submission['outputs_dir'], outputs_key)
        shutil.move(src, outputs_dest)

def run_job_submission(*args, JobSubmissionRunner=None, job=None, cfg=None,
                       submission=None, **kwargs):
    return JobSubmissionRunner(*args, job=job, cfg=cfg, submission=submission,
                               **kwargs).run_submission()
