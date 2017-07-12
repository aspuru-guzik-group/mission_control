import os
import subprocess

from mc.job_module_utils.jobdir_runner import JobdirRunner

from .workdir_builder import WorkdirBuilder


class ExampleEchoJobdirRunner(JobdirRunner):
    class InvalidMessageError(Exception): pass

    def _run_jobdir(self):
        self._validate_job_params()
        workdir_meta = self._build_workdir()
        self._run_workdir(workdir_meta=workdir_meta)
        workdir_output = os.path.join(workdir_meta['dir'],
                                      workdir_meta['output_file_name'])
        self.move_to_outputs(src=workdir_output, outputs_key='example_echo.out')

    def _validate_job_params(self):
        try: assert self.job_params['message'] is not None
        except Exception as exc: raise self.InvalidMessageError()

    def _build_workdir(self):
        workdir_meta = WorkdirBuilder(
            workdir=self.mkdtemp(prefix='example_echo.wd.'),
            workdir_params=self.get_workdir_params()
        ).build_workdir()
        return workdir_meta

    def get_workdir_params(self):
        workdir_params = {'message': self.job_params.get('message')}
        return workdir_params

    def _run_workdir(self, workdir_meta=None):
        entrypoint_path = os.path.join(workdir_meta['dir'],
                                       workdir_meta['entrypoint'])
        cmd = [entrypoint_path]
        extra_env_vars = self.cfg.get('example_echo', {}).get('env_vars', {})
        workdir_env = {**os.environ, **extra_env_vars}
        subprocess.run(cmd, env=workdir_env, check=True)

run_jobdir = ExampleEchoJobdirRunner.run_jobdir

