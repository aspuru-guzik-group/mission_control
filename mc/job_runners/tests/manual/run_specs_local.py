import json
import os
import sys
import tempfile
import textwrap
from uuid import uuid4

from mc.local_job_runners.base_daemon import BaseDaemon
from mc.local_job_runners.job_client import JobClient


def main():
    daemon = BaseDaemon(job_client=generate_job_client(),
                        job_dir_factory=generate_job_dir_factory(),
                        tick_interval=3)
    daemon.start()

def generate_job_client():
    job_paths = sys.argv[1:]
    jobs = load_jobs_from_paths(job_paths)
    job_client = SimpleJobClient(jobs=jobs)
    return job_client

def load_jobs_from_paths(job_paths):
    jobs = []
    for job_path in job_paths:
        job = json.load(open(job_path))
        job.setdefault('uuid', str(uuid4()))
        jobs.append(job)
    return jobs

class SimpleJobClient(JobClient):
    def __init__(self, jobs=None):
        self.jobs = {job['uuid']: job for job in jobs}
        self.claimed_uuids = {}

    def fetch_jobs(self):
        claimable_jobs = []
        for uuid, job in self.jobs.items():
            if uuid not in self.claimed_uuids and \
               job.get('status') is not self.Statuses.COMPLETED:
                claimable_jobs.append(job)
        return claimable_jobs

    def claim_jobs(self, uuids=None):
        for _uuid in uuids:
            self.claimed_uuids[_uuid] = True
        return {_uuid: True for _uuid in uuids}

    def update_jobs(self, updates_by_uuid=None):
        for uuid, updates_for_uuid in updates_by_uuid.items():
            self.jobs[uuid].update(updates_for_uuid)

def generate_job_dir_factory():
    return StubJobDirFactory()

class StubJobDirFactory():
    def build_dir_for_spec(self, job=None):
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

if __name__ == '__main__':
    main()
