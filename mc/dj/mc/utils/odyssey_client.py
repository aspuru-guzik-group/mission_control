from datetime import datetime
import os
import re
import sys
import time

REMOTE_JOBS_DIR = '~/.odyssey_job_runner/jobs'

class OdysseyClient(object):
    def __init__(self, ssh_client=None):
        self.ssh_client = ssh_client

    def run_job(self, local_job_dir):
        job_metadata = self.start_job(local_job_dir)
        job_state = self.get_job_state(job_metadata)
        while job_state['JobState'] in ['PENDING', 'RUNNING']:
            time.sleep(30)
            job_state = self.get_job_state(job_metadata)
        print("Job completed.", file=sys.stderr)

    def start_job(self, local_job_dir):
        job_metadata = {'local_dir': local_job_dir}
        job_metadata['remote_dir'] = self._upload_job_dir(local_job_dir)
        job_cmd = os.path.join(job_metadata['remote_dir'], 'job.sh')
        odyssey_cmd = ['sbatch', job_cmd]
        submission_output = self.ssh_client.run(
            odyssey_cmd, check=True).stdout
        match = re.match(r'Submitted batch job (\d+)', submission_output)
        if match:
            job_metadata['job_id'] = match.group(1)
        return job_metadata

    def _upload_job_dir(self, local_job_dir):
        basename = os.path.basename(local_job_dir)
        self._ensure_remote_jobs_dir()
        remote_job_dir = os.path.join(
            REMOTE_JOBS_DIR, '%s.%s' % (basename, datetime.now().timestamp()))
        self.ssh_client.scp(
            src=local_job_dir,
            dst=remote_job_dir,
            flags='-r')
        return remote_job_dir

    def _ensure_remote_jobs_dir(self):
        remote_jobs_dir_exists = self.ssh_client.run(
            ['[ -e %s ]' % REMOTE_JOBS_DIR], check=False).returncode == 0
        if not remote_jobs_dir_exists:
            self.ssh_client.run(['mkdir -p %s' % REMOTE_JOBS_DIR])

    def get_job_state(self, job_metadata):
        odyssey_cmd = ['scontrol', 'show --details --oneliner',
                       'jobid %s' % job_metadata['job_id']]
        raw_job_state = self.ssh_client.run(odyssey_cmd).stdout
        return self._parse_raw_job_state(raw_job_state)

    def _parse_raw_job_state(self, raw_job_state):
        job_state = {}
        for key_value_str in raw_job_state.split():
            key, value = key_value_str.split('=', 1)
            job_state[key] = value
        return job_state
