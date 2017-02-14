import os
from .slurm_execution_client import SlurmExecutionClient

class RemoteSlurmExecutionClient(object):
    def __init__(self, ssh_client=None, remote_workdir=None):
        self.ssh_client = ssh_client
        self.remote_workdir = '.remote_slurm_execution_client'
        self.slurm_execution_client = SlurmExecutionClient(
            process_runner=self.ssh_client)

    def start_execution(self, job=None):
        remote_dir_meta = self.upload_job(job=job)
        # Execute job, but use remote dir instead of local dir.
        execution_meta = self.slurm_execution_client.start_execution(
            job={**job, 'dir': {**job['dir'], 'dir': remote_dir_meta['dir']}})
        execution_meta['remote_dir'] = remote_dir_meta
        return execution_meta

    def upload_job(self, job=None):
        self.ensure_remote_workdir()
        local_src = job['dir']['dir']
        remote_dest = os.path.join(self.remote_workdir, job['uuid'])
        self.ssh_client.rsync_to_remote(local_src_path=local_src + '/',
                                        remote_dest_path=remote_dest,
                                        flags='-a')
        remote_dir_meta = {'dir': remote_dest}
        return remote_dir_meta

    def ensure_remote_workdir(self):
        self.ssh_client.run_process(cmd=['mkdir', '-p', self.remote_workdir],
                                    check=True)

    def get_execution_state(self, job=None):
        execution_state = self.slurm_execution_client.get_execution_state(
            job=job)
        if not execution_state['executing']:
            self.on_job_completed(job=job)
        return execution_state

    def on_job_completed(self, job=None):
        self.download_job(job=job)

    def download_job(self, job=None):
        remote_src = job['execution']['remote_dir']
        local_dest = job['dir']['dir']
        self.ssh_client.rsync_from_remote(remote_src_path=remote_src + '/',
                                          local_dest_path=local_dest,
                                          flags='-a')
