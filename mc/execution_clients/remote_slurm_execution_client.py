import os
from uuid import uuid4

from .slurm_execution_client import SlurmExecutionClient



class RemoteSlurmExecutionClient(object):
    def __init__(self, ssh_client=None, remote_workdir=None):
        self.ssh_client = ssh_client
        self.rsync_flags = ['-a', '--no-owner']
        self.remote_workdir = '.remote_slurm_execution_client'
        self.slurm_execution_client = SlurmExecutionClient(
            process_runner=self.ssh_client)

    def start_execution(self, submission=None):
        remote_dir_meta = self.upload_submission(submission=submission)
        submission_w_remote_dir = {**submission, 'dir': remote_dir_meta['dir']}
        execution_meta = self.slurm_execution_client.start_execution(
            submission=submission_w_remote_dir)
        execution_meta['remote_dir'] = remote_dir_meta
        return execution_meta

    def upload_submission(self, submission=None):
        self.ensure_remote_workdir()
        remote_dir_path = self.generate_remote_dir_path()
        self.ssh_client.rsync_to_remote(
            local_src_path=(submission['dir'] + '/'),
            remote_dest_path=remote_dir_path,
            flags=self.rsync_flags
        )
        remote_dir_meta = {'dir': remote_dir_path}
        return remote_dir_meta

    def generate_remote_dir_path(self):
        return os.path.join(self.remote_workdir, str(uuid4()))

    def ensure_remote_workdir(self):
        self.ssh_client.run_process(cmd=['mkdir', '-p', self.remote_workdir],
                                    check=True)

    def get_execution_state(self, execution_meta=None):
        execution_state = self.slurm_execution_client.get_execution_state(
            execution_meta=execution_meta)
        if execution_state['run_status'] != 'RUNNING':
            self.on_execution_completed(execution_meta=execution_meta)
        return execution_state

    def on_execution_completed(self, execution_meta=None):
        self.download_completed_dir(execution_meta=execution_meta)

    def download_completed_dir(self, execution_meta=None):
        self.ssh_client.rsync_from_remote(
            remote_src_path=execution_meta['remote_dir']['dir'] + '/',
            local_dest_path=execution_meta['submission']['dir'],
            flags=self.rsync_flags
        )
