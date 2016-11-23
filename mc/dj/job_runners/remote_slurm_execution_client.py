import os
import shutil
from .slurm_execution_client import SlurmExecutionClient

class RemoteSlurmExecutionClient(object):
    def __init__(self, ssh_client=None, remote_workdir=None):
        self.ssh_client = ssh_client
        self.remote_workdir = '.remote_slurm_execution_client'
        self.slurm_execution_client = SlurmExecutionClient(
            process_runner=self.ssh_client)
        self.remote_dir_key = '_%s:remote_dir' % id(self)

    def start_execution(self, job=None):
        remote_dir_meta = self.upload_job(job=job)
        # Execute job, but use remote dir instead of local dir.
        executing_job = self.slurm_execution_client.start_execution(
            job={**job, 'dir': {**job['dir'], 'dir': remote_dir_meta['dir']}})
        # Return modified executing job, so that it has its original 'dir'
        # value.
        wrapped_executing_job = {**executing_job, 'dir': job['dir'],
                                 self.remote_dir_key: remote_dir_meta}
        return wrapped_executing_job

    def upload_job(self, job=None):
        self.ensure_remote_workdir()
        local_src = job['dir']['dir']
        remote_dest = os.path.join(self.remote_workdir, job['uuid'])
        self.ssh_client.scp_to(src=local_src, dest=remote_dest, flags='-r')
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
        self.rm_orig_job_dir(job=job)
        self.download_job(job=job)

    def rm_orig_job_dir(self, job=None):
        shutil.rmtree(job['dir']['dir'])

    def download_job(self, job=None):
        remote_src = job['remote_dir']['dir']
        local_dest = job['dir']['dir']
        self.ssh_client.scp_from(src=remote_src, dest=local_dest, flags='-r')
