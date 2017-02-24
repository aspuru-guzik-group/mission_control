import os

from mc.job_runners.base_job_runner import BaseJobRunner
from mc.job_runners.remote_slurm_execution_client \
        import RemoteSlurmExecutionClient as ExecutionClient

class OdysseyPushJobRunner:
    def __init__(self, *args, setup=True, **base_job_runner_kwargs):
        if setup: self.setup(*args, **base_job_runner_kwargs)

    def setup(self, *args, ssh_client=None, execution_client=None,
              **base_job_runner_kwargs):
        execution_client = execution_client or self.generate_execution_client(
            ssh_client=ssh_client)
        self.base_job_runner = BaseJobRunner(
            execution_client=execution_client,
            get_job_execution_result=self.get_job_execution_result,
            **base_job_runner_kwargs
        )
    
    def generate_execution_client(self, ssh_client=None):
        return ExecutionClient(ssh_client=ssh_client)

    def get_job_execution_result(self, job_state=None):
        def get_checkpoint_path(checkpoint_name=None):
            return os.path.join(
                job_state['completed_dir'],
                job_state['submission']['checkpoint_files']\
                [checkpoint_name]
            )
        if os.path.exists(get_checkpoint_path(checkpoint_name='completed')):
            result = {'result': 'COMPLETED'}
        else:
            failure_checkpoint_path = get_checkpoint_path(
                checkpoint_name='failed')
            try:
                error = open(failure_checkpoint_path).read()
            except Exception:
                error = ("Cause unknown, could not read failure checkpoint file" 
                         " '{}'.").format(failure_checkpoint_path)
            result = {'result': 'FAILED', 'error': error}
        return result

    def run(self, *args, **kwargs):
        self.base_job_runner.run(*args, **kwargs)

    def tick(self, *args, **kwargs):
        self.base_job_runner.tick(*args, **kwargs)
