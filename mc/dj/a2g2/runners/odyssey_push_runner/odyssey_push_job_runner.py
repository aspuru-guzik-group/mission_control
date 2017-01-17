from job_runners.base_job_runner import BaseJobRunner
from job_runners.remote_slurm_execution_client import RemoteSlurmExecutionClient
from job_runners.ssh_control_socket_client import SSHControlSocketClient


class OdysseyPushJobRunner:
    def __init__(self, *args, setup=True, odyssey_user=None, odyssey_host=None,
                 **base_job_runner_kwargs):
        self.odyssey_user = odyssey_user
        self.odyssey_host = odyssey_host
        if setup: self.setup(*args, **base_job_runner_kwargs)

    def setup(self, *args, execution_client=None,
              transfer_client=None, **base_job_runner_kwargs):
        execution_client = execution_client or self.generate_execution_client()
        transfer_client = transfer_client or self.generate_transfer_client()
        self.base_job_runner = BaseJobRunner(
            execution_client=execution_client,
            transfer_client=transfer_client,
            **base_job_runner_kwargs
        )
    
    def generate_execution_client(self):
        ssh_client = SSHControlSocketClient(user=self.odyssey_user,
                                            host=self.odyssey_host)
        return RemoteSlurmExecutionClient(ssh_client=ssh_client)

    def generate_transfer_client(self):
        class StubTransferClient: pass
        return StubTransferClient()

    def run(self, *args, **kwargs):
        self.base_job_runner.run(*args, **kwargs)

    def tick(self, *args, **kwargs):
        self.base_job_runner.tick(*args, **kwargs)
