import psutil
import subprocess
import time


class LocalJobRunnerDaemon(object):
    def __init__(self, job_client=None, job_dir_factory=None, interval=120,
                 max_running_jobs=3):
        self.job_client = job_client
        self.job_dir_factory = job_dir_factory
        self.interval = interval
        self.max_running_jobs = max_running_jobs

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
        partial_job = {'job_spec': job_spec, 'dir': dir_meta}
        spec_meta = self.run_job(job=partial_job)
        full_job = {**partial_job, 'proc': spec_meta}
        self.running_jobs[job_spec['uuid']] = full_job

    def claim_job_spec(self, job_spec=None):
        claim_results = self.job_client.claim_jobs(uuids=[job_spec['uuid']])
        return claim_results.get(job_spec['uuid'], False)

    def build_job_dir(self, job_spec=None):
        job_dir_meta = self.job_dir_factory.build_job_dir(job_spec=job_spec)
        return job_dir_meta

    def run_job(self, job=None):
        cmd = 'cd {dir}; {entrypoint}'.format(**job['dir'])
        proc = subprocess.Popen(cmd, shell=True)
        return {'pid': proc.pid}

    def process_running_jobs(self):
        current_running_job_states = self.get_running_job_states()
        completed_jobs = [job for job_uuid, job in self.running_jobs.items()
                          if current_running_job_states[job_uuid] is None]
        for completed_job in completed_jobs:
            self.process_completed_job(job=completed_job)

    def get_running_job_states(self):
        job_states = {}
        for job_uuid, running_job in self.running_jobs.items():
            try:
                job_state = psutil.Process(
                    pid=running_job['proc']['pid']).as_dict()
            except psutil.NoSuchProcess:
                job_state = None
            job_states[job_uuid] = job_state
        return job_states

    def process_completed_job(self, job=None):
        transfer_meta = self.transfer_job(job=job)
        self.transferring_jobs[job['job_spec']['uuid']] = {
            **job, 'transfer': transfer_meta}
        del self.running_jobs[job['job_spec']['uuid']]

    def transfer_job(self, job=None, job_dir_meta=None):
        pass

    def process_transferring_jobs(self):
        pass
