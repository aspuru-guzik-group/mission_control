import json
import os
import sys
import tempfile
import textwrap
from uuid import uuid4
from .base_daemon import BaseDaemon


def main():
    daemon = BaseDaemon(job_spec_client=generate_job_spec_client(),
                        job_dir_factory=generate_job_dir_factory(),
                        transfer_client=generate_transfer_client())
    daemon.run()

def generate_job_spec_client():
    job_spec_paths = sys.argv[1:]
    job_specs = load_job_specs_from_paths(job_spec_paths)
    job_spec_client = SimpleJobSpecClient(job_specs=job_specs)
    return job_spec_client

def load_job_specs_from_paths(job_spec_paths):
    job_specs = {}
    for job_spec_path in job_spec_paths:
        job_spec = json.load(open(job_spec_path))
        job_spec.setdefault('uuid', uuid4())
    return job_specs

class SimpleJobSpecClient(object):
    def __init__(self, job_specs=None):
        self.job_specs = {job_spec['uuid']: job_spec for job_spec in job_specs}
        self.claimed_uuids = {}

    def fetch_job_specs(self):
        claimable_job_specs = []
        for uuid, job_spec in self.job_specs:
            if uuid not in self.claimed_uuids and \
               job_spec['status'] is not self.statuses.COMPLETED:
                claimable_job_specs.append(job_spec)
        return claimable_job_specs

    def claim_job_specs(self, uuids=None):
        for uuid in uuids:
            self.claimed_uuids[uuid] = True

    def update_job_specs(self, updates_by_uuid=None):
        for uuid, updates_for_uuid in updates_by_uuid.items():
            self.job_specs[uuid].update(updates_for_uuid)

def generate_job_dir_factory():
    return StubJobDirFactory()

class StubJobDirFactory():
    def build_dir_for_spec(job_spec=None):
        dir_path = tempfile.mkdtemp(prefix='stjdf.')
        entrypoint_path = os.path.join(dir_path, 'job.sh')
        with open(entrypoint_path, 'w') as f:
            f.write(textwrap.dedent(
                """
                #!/bin/bash
                echo beef > meat.txt
                """))
        dir_meta = {
            'dir': dir_path,
            'entrypoint': entrypoint_path
        }
        return dir_meta

def generate_transfer_client():
    return StubTransferClient()

class StubTransferClient():
    def start_transfer(self, job=None):
        storage_uri = job['job_spec']['storage']['uri']
        storage_client = self.get_storage_client_for_storage_uri(storage_uri)
        key = storage_client.put(job['dir']['dir'])
        return {'key': key}

    def get_storage_client_for_storage_uri(self, storage_uri):
        pass

    def get_transfer_state(self, key=None):
        pass

if __name__ == '__main__':
    main()
