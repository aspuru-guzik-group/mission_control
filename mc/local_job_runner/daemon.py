import time

class LocalJobRunnerDaemon(object):
    def __init__(self, job_client=None, interval=120, max_concurrent_jobs=3):
        self.job_client = job_client
        self.interval = interval
        self.max_concurrent_jobs = max_concurrent_jobs

        self._ticking = False
        self.running_jobs = {}

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
        num_job_slots = self.max_concurrent_jobs - len(self.running_jobs)
        if num_job_slots <= 0: return
        candidate_jobs = self.fetch_candidate_jobs()
        for job in candidate_jobs[:-num_job_slots]:
            self.process_job(job)
        self.poll_running_jobs()

    def fetch_candidate_jobs(self):
        return self.job_client.fetch_jobs()

    def process_job(self, job):
        claimed = self.claim_job(job)
        if not claimed: return
        job_dir = self.build_job_dir(job=job)
        job_proc_meta = self.run_job(job=job, job_dir=job_dir)
        self.running_jobs[job['uuid']] = {'job': job, 'run_dir': job_dir,
                                          'proc': job_proc_meta}

    def claim_job(self, job=None):
        claim_results = self.job_client.claim_jobs(uuids=[job['uuid']])
        return claim_results.get(job['uuid'], False)

    def poll_running_jobs(self):
        pass

    def build_job_dir(self):
        pass

    def run_job(self, job_meta=None):
        pass

    def get_running_jobs(self):
        return self.running_jobs
