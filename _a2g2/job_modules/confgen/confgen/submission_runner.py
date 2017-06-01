import os
import subprocess

from ...a2g2_common.base_submission_runner import BaseSubmissionRunner
from .workdir_builder import WorkdirBuilder
from .. import constants as confgen_constants


class SubmissionRunner(BaseSubmissionRunner):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run_submission(self):
        workdir_meta = self.create_workdir()
        try:
            self.run_workdir(workdir_meta=workdir_meta)
        finally:
            self.move_to_outputs(
                src=workdir_meta['dir'],
                outputs_key=confgen_constants.CONFGEN_OUTPUTS_KEY)

    def create_workdir(self):
        workdir_builder = WorkdirBuilder(
            workdir=self.generate_tmp_dir(),
            workdir_params=self.get_workdir_params())
        workdir_meta = workdir_builder.build_workdir()
        return workdir_meta

    def get_workdir_params(self):
        workdir_params = {
            'confgen_params': self.job['job_spec']['job_params'].get(
                'confgen_params')
        }
        return workdir_params

    def run_workdir(self, workdir_meta=None):
        entrypoint_path = os.path.join(workdir_meta['dir'],
                                       workdir_meta['entrypoint'])
        cmd = ['/bin/bash', entrypoint_path]
        workdir_env = {
            **os.environ, 
            **(self.cfg.get('a2g2.jobs.confgen.confgen', {}).get(
                'env_vars', {}))
        }
        subprocess.run(cmd, env=workdir_env, check=True)

def run_job_submission(*args, job=None, cfg=None, submission=None, **kwargs):
    runner = SubmissionRunner(job=job, cfg=cfg, submission=submission)
    return runner.run_submission()
