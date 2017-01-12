import os
import tempfile
import textwrap
from uuid import uuid4

from job_runners.base_job_runner import BaseJobRunner
from job_runners.filesystem_storage_client import FileSystemStorageClient
from job_runners.job_client import JobClient
from job_runners.remote_slurm_execution_client import RemoteSlurmExecutionClient
from job_runners.ssh_control_socket_client import SSHControlSocketClient


def main():
    runner = BaseJobRunner(
        job_client=generate_job_client(),
        job_dir_factory=generate_job_dir_factory(),
        execution_client=generate_execution_client(),
        transfer_client=generate_transfer_client(),
        tick_interval=3
    )
    runner.start()

def generate_job_client():
    jobs = [{'uuid': str(uuid4())} for i in range(3)]
    job_client = StubJobClient(jobs=jobs)
    return job_client

class StubJobClient(JobClient):
    def __init__(self, jobs=None):
        self.jobs = {job['uuid']: job for job in jobs}
        self.claimed_uuids = {}

    def fetch_claimable_jobs(self):
        claimable_jobs = []
        for uuid, job in self.jobs.items():
            if uuid not in self.claimed_uuids and \
               job.get('status') is not self.Statuses.Completed:
                claimable_jobs.append(job)
        return claimable_jobs

    def claim_jobs(self, uuids=None):
        claimed_specs = {}
        for _uuid in uuids:
            self.claimed_uuids[_uuid] = True
            claimed_specs[_uuid] = self.jobs[_uuid]
        return claimed_specs

    def update_jobs(self, updates_by_uuid=None):
        for uuid, updates_for_uuid in updates_by_uuid.items():
            self.jobs[uuid].update(updates_for_uuid)

def generate_job_dir_factory():
    return StubJobDirFactory()

class StubJobDirFactory():
    def build_dir_for_spec(self, job=None):
        dir_path = tempfile.mkdtemp(prefix='stjdf.')
        entrypoint_file_name = 'job.sh'
        entrypoint_path = os.path.join(dir_path, entrypoint_file_name)
        with open(entrypoint_path, 'w') as f:
            f.write(textwrap.dedent("""
                #!/bin/bash
                echo beef > meat.txt
                """).strip())
        os.chmod(entrypoint_path, 0o775)
        dir_meta = {
            'dir': dir_path,
            'entrypoint': entrypoint_file_name,
        }
        return dir_meta

def generate_execution_client():
    ssh_client = SSHControlSocketClient(user=os.getenv('SLURM_USER'),
                                        host=os.getenv('SLURM_HOST'))
    ssh_client.connect()
    execution_client = RemoteSlurmExecutionClient(ssh_client=ssh_client)
    return execution_client

def generate_transfer_client():
    return StubTransferClient()

class StubTransferClient(object):
    def __init__(self):
        self.state_dir = tempfile.mkdtemp(prefix='stc.state')

    def start_transfer(self, job=None):
        storage_client = self.get_storage_client_for_job(job=job)
        storage_key = storage_client.put(job['dir']['dir'])
        transfer_meta = {'key': storage_key}
        return transfer_meta

    def get_storage_client_for_job(self, job=None):
        storage_client = FileSystemStorageClient(
            storage_root_path='/tmp/TEST_STORAGE',
            state_file_path=os.path.join(self.state_dir,
                                         'TEST_STORAGE.json')
        )
        return storage_client

    def get_transfer_state(self, job=None):
        storage_client = self.get_storage_client_for_job(job=job)
        return storage_client.poll(key=job['transfer']['key'])

if __name__ == '__main__':
    main()
