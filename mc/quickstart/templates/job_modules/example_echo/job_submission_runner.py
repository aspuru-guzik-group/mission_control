import os
import subprocess

from mc.job_module_utils.base_job_submission_runner \
        import BaseJobSubmissionRunner

from .workdir_builder import WorkdirBuilder


class JobSubmissionRunner(BaseJobSubmissionRunner):
    class InvalidMessageError(Exception): pass

    def validate_job_params(self, job_params=None):
        job_params = job_params or self.job_params
        try: assert job_params['message'] is not None
        except Exception as exc: raise self.InvalidMessageError()

    def _run_job_submission(self):
        workdir_meta = self.create_workdir()
        self.run_workdir(workdir_meta=workdir_meta)
        self.move_to_outputs(
            src=os.path.join(workdir_meta['dir'],
                             workdir_meta['output_file_name']),
            outputs_key='example_echo.output'
        )

    def create_workdir(self):
        workdir_meta = WorkdirBuilder(
            workdir=self.generate_tmp_dir(prefix='example_echo.wd.'),
            workdir_params=self.get_workdir_params()
        ).build_workdir()
        return workdir_meta

    def get_workdir_params(self):
        workdir_params = {'message': self.job_params.get('message')}
        return workdir_params

    def run_workdir(self, workdir_meta=None):
        entrypoint_path = os.path.join(workdir_meta['dir'],
                                       workdir_meta['entrypoint'])
        cmd = [entrypoint_path]
        extra_env_vars = self.cfg.get('example_echo', {}).get('env_vars', {})
        workdir_env = {**os.environ, **extra_env_vars}
        subprocess.run(cmd, env=workdir_env, check=True)
