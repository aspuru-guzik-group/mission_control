import os
import tempfile
import time
import shutil
import subprocess

from .workdir_builder import WorkdirBuilder


class SubmissionRunner(object):
    def __init__(self, job=None, cfg=None, submission=None):
        self.job = job
        self.cfg = cfg
        self.submission = submission

    def run_submission(self):
        workdir_meta = self.create_workdir()
        try:
            self.run_workdir(workdir_meta=workdir_meta)
        finally:
            self.move_workdir_to_outputs(workdir_meta=workdir_meta)

    def create_workdir(self):
        workdir_builder = WorkdirBuilder(
            workdir=self.generate_workdir(),
            workdir_params=self.get_workdir_params())
        workdir_meta = workdir_builder.build_workdir()
        return workdir_meta

    def generate_workdir(self):
        scratch_dir = self.submission.get('scratch_dir') or tempfile.mkdtemp()
        workdir = os.path.join(scratch_dir, 'confgen.%s' % time.time())
        os.makedirs(workdir)
        return workdir

    def get_workdir_params(self):
        workdir_params = {
            'confgen_params': self.job['job_spec']['job_params'].get(
                'confgen_params')
        }
        return workdir_params

    def run_workdir(self, workdir_meta=None):
        cmd = ['bash', workdir_meta['entrypoint']]
        workdir_env = {
            **os.environ,
            'CONFGEN_EXE': self.submission.get('ctx', {}).get('CONFGEN_EXE')
        }
        subprocess.run(cmd, check=True, env=workdir_env)

    def move_workdir_to_outputs(self, workdir_meta=None):
        outputs_dest = os.path.join(self.submission['dir']['outputs'],
                                    'confgen')
        shutil.move(workdir_meta['dir'], outputs_dest)
