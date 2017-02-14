from job_runners.base_job_runner import BaseJobRunner
from job_runners.remote_slurm_execution_client import RemoteSlurmExecutionClient


class OdysseyPushJobRunner:
    def __init__(self, *args, setup=True, **base_job_runner_kwargs):
        if setup: self.setup(*args, **base_job_runner_kwargs)

    def setup(self, *args, ssh_client=None, execution_client=None,
              **base_job_runner_kwargs):
        execution_client = execution_client or self.generate_execution_client(
            ssh_client=ssh_client)
        self.base_job_runner = BaseJobRunner(execution_client=execution_client,
                                             **base_job_runner_kwargs)
    
    def generate_execution_client(self, ssh_client=None):
        return RemoteSlurmExecutionClient(ssh_client=ssh_client)

    def run(self, *args, **kwargs):
        self.base_job_runner.run(*args, **kwargs)

    def tick(self, *args, **kwargs):
        self.base_job_runner.tick(*args, **kwargs)
