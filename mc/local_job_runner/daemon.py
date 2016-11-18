import psutil
import subprocess
import time


class LocalJobRunnerDaemon(object):
    def __init__(self, job_client=None, job_dir_factory=None, interval=120,
                 max_running_jobs=3, transfer_client=None):
        self.job_client = job_client
        self.job_dir_factory = job_dir_factory
        self.interval = interval
        self.max_running_jobs = max_running_jobs
        self.transfer_client = transfer_client

        self._ticking = False
        self.running_jobs = {}
        self.transferring_jobs = {}

    def start(self):
        self._ticking = True
        self.run()

    def stop(self):
        self._ticking = False

    def run(self, ntimes=None):
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep()
        else:
            while self._ticking:
                self._tick_and_sleep()

    def _tick_and_sleep(self):
        self.tick()
        time.sleep(self.interval)

    def tick(self):
        num_job_slots = self.max_running_jobs - len(self.running_jobs)
        if num_job_slots <= 0: return
        candidate_job_specs = self.fetch_candidate_job_specs()
        for candidate_job_spec in candidate_job_specs[:-num_job_slots]:
            self.process_candidate_job_spec(job_spec=candidate_job_spec)
        self.process_running_jobs()
        self.process_transferring_jobs()

    def fetch_candidate_job_specs(self):
        return self.job_client.fetch_jobs()

    def process_candidate_job_spec(self, job_spec=None):
        claimed = self.claim_job_spec(job_spec)
        if not claimed: return
        dir_meta = self.build_job_dir(job_spec=job_spec)
        partial_job = {'key': job_spec['uuid'], 'job_spec': job_spec,
                       'dir': dir_meta}
        spec_meta = self.start_job_run(job=partial_job)
        full_job = {**partial_job, 'proc': spec_meta}
        self.running_jobs[partial_job['key']] = full_job

    def claim_job_spec(self, job_spec=None):
        claim_results = self.job_client.claim_jobs(uuids=[job_spec['uuid']])
        return claim_results.get(job_spec['uuid'], False)

    def build_job_dir(self, job_spec=None):
        job_dir_meta = self.job_dir_factory.build_job_dir(job_spec=job_spec)
        return job_dir_meta

    def start_job_run(self, job=None):
        cmd = 'cd {dir}; {entrypoint}'.format(**job['dir'])
        proc = subprocess.Popen(cmd, shell=True)
        return {'pid': proc.pid}

    def process_running_jobs(self):
        job_run_states = self.get_job_run_states()
        completed_jobs = [job for job_key, job in self.running_jobs.items()
                          if job_run_states[job_key] is None]
        for completed_job in completed_jobs:
            self.process_completed_job(job=completed_job)

    def get_job_run_states(self):
        job_run_states = {}
        for job_key, job in self.running_jobs.items():
            try:
                job_run_state = psutil.Process(
                    pid=job['proc']['pid']).as_dict()
            except psutil.NoSuchProcess:
                job_run_state = None
            job_run_states[job_key] = job_run_state
        return job_run_states

    def process_completed_job(self, job=None):
        transfer_meta = self.start_job_transfer(job=job)
        self.transferring_jobs[job['key']] = {**job, 'transfer': transfer_meta}
        del self.running_jobs[job['key']]

    def start_job_transfer(self, job=None):
        transfer_meta = self.transfer_client.start_transfer(job=job)
        return transfer_meta

    def process_transferring_jobs(self):
        job_transfer_states = self.get_job_transfer_states()
        transferred_jobs = [
            job for job_key, job in self.transferring_jobs.items()
            if job_transfer_states[job_key] is None
        ]
        for transferred_job in transferred_jobs:
            self.process_transferred_job(job=transferred_job)
        pass

    def get_job_transfer_states(self):
        job_transfer_states = {}
        for job_key, job in self.transferring_jobs.items():
            job_transfer_state = self.transfer_client.get_transfer_state(
                key=job['transfer']['key'])
            job_transfer_states[job_key] = job_transfer_state
        return job_transfer_states

    def process_transferred_job(self, job=None):
        self.update_job_spec(job_spec=job['job_spec'], updates={
            'status': self.job_client.statuses.TRANSFERRED,
            'transfer_meta': job['transfer'],
        })
        del self.transferring_jobs[job['key']]

    def update_job_spec(self, job_spec=None, updates=None):
        self.job_client.update_job(uuid=job_spec['uuid'], updates=updates)

