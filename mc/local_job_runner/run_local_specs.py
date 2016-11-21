import json
import os
import sys
import tempfile
import textwrap
from uuid import uuid4
from .base_daemon import BaseDaemon
from .filesystem_storage_client import FileSystemStorageClient
from .job_spec_client import JobSpecClient


def main():
    daemon = BaseDaemon(job_spec_client=generate_job_spec_client(),
                        job_dir_factory=generate_job_dir_factory(),
                        transfer_client=generate_transfer_client(),
                        tick_interval=3)
    daemon.start()

def generate_job_spec_client():
    job_spec_paths = sys.argv[1:]
    job_specs = load_job_specs_from_paths(job_spec_paths)
    job_spec_client = SimpleJobSpecClient(job_specs=job_specs)
    return job_spec_client

def load_job_specs_from_paths(job_spec_paths):
    job_specs = []
    for job_spec_path in job_spec_paths:
        job_spec = json.load(open(job_spec_path))
        job_spec.setdefault('uuid', str(uuid4()))
        job_specs.append(job_spec)
    return job_specs

class SimpleJobSpecClient(JobSpecClient):
    def __init__(self, job_specs=None):
        self.job_specs = {job_spec['uuid']: job_spec for job_spec in job_specs}
        self.claimed_uuids = {}

    def fetch_job_specs(self):
        claimable_job_specs = []
        for uuid, job_spec in self.job_specs.items():
            if uuid not in self.claimed_uuids and \
               job_spec.get('status') is not self.Statuses.COMPLETED:
                claimable_job_specs.append(job_spec)
        return claimable_job_specs

    def claim_job_specs(self, uuids=None):
        for _uuid in uuids:
            self.claimed_uuids[_uuid] = True
        return {_uuid: True for _uuid in uuids}

    def update_job_specs(self, updates_by_uuid=None):
        for uuid, updates_for_uuid in updates_by_uuid.items():
            self.job_specs[uuid].update(updates_for_uuid)

def generate_job_dir_factory():
    return StubJobDirFactory()

class StubJobDirFactory():
    def build_dir_for_spec(self, job_spec=None):
        dir_path = tempfile.mkdtemp(prefix='stjdf.')
        entrypoint_path = os.path.join(dir_path, 'job.sh')
        with open(entrypoint_path, 'w') as f:
            f.write(textwrap.dedent(
                """
                #!/bin/bash
                echo beef > meat.txt
                """))
        os.chmod(entrypoint_path, 0o775)
        dir_meta = {
            'dir': dir_path,
            'entrypoint': entrypoint_path
        }
        return dir_meta

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
        storage_uri = job['job_spec']['storage']['uri']
        return self.get_storage_client_for_storage_uri(storage_uri)

    def get_storage_client_for_storage_uri(self, storage_uri):
        if storage_uri == 'TEST_STORAGE':
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
