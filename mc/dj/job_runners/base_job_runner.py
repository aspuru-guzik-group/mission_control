import logging
import time


class BaseJobRunner(object):
    def __init__(self, job_client=None, job_dir_factory=None,
                 execution_client=None, tick_interval=120, max_executing_jobs=3,
                 logger=None):
        self.job_client = job_client
        self.job_dir_factory = job_dir_factory
        self.execution_client = execution_client
        self.tick_interval = tick_interval
        self.max_executing_jobs = max_executing_jobs
        self.logger = logger or logging

        self._ticking = False
        self.tick_counter = 0
        self.jobs = {}
        self.executing_jobs = {}

    def start(self):
        self._ticking = True
        self.run()

    def stop(self):
        self._ticking = False

    def run(self, ntimes=None, tick_interval=None):
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep(tick_interval=tick_interval)
        else:
            while self._ticking:
                self._tick_and_sleep(tick_interval=tick_interval)

    def _tick_and_sleep(self, tick_interval=None):
        if tick_interval is None:
            tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        self.tick_counter += 1
        logging.debug('tick #%s' % self.tick_counter)
        self.process_executing_jobs()
        num_job_slots = self.max_executing_jobs - len(self.executing_jobs)
        if num_job_slots <= 0: return
        claimable_jobs = self.fetch_claimable_jobs()
        for claimable_job in claimable_jobs[:num_job_slots]:
            self.process_claimable_job(job=claimable_job)

    def fetch_claimable_jobs(self):
        logging.debug('fetch_claimable_jobs')
        return self.job_client.fetch_claimable_jobs()

    def process_claimable_job(self, job=None):
        logging.debug('process_claimable_job')
        claimed_job = self.claim_job(job)
        if not claimed_job: return
        self.jobs[claimed_job['uuid']] = claimed_job
        try:
            claimed_job['dir'] = self.build_job_dir(job=claimed_job)
            claimed_job['execution'] = self.start_job_execution(job=claimed_job)
        except Exception as error:
            self.fail_job(job=claimed_job, error=error)

    def fail_job(self, job=None, error=None):
        logging.exception(error)
        self.update_job(job=job, updates={'status': 'FAILED', 'error': error})
        del self.jobs[job['uuid']]

    def claim_job(self, job=None):
        logging.debug('claim_job')
        claimed_jobs = self.job_client.claim_jobs(
            uuids=[job['uuid']])
        return claimed_jobs.get(job['uuid'], False)

    def build_job_dir(self, job=None):
        logging.debug('build_job_dir')
        job_dir_meta = self.job_dir_factory.build_dir_for_job(
            job=job)
        return job_dir_meta

    def start_job_execution(self, job=None):
        logging.debug('start_job_execution')
        self.executing_jobs[job['uuid']] = job
        execution_meta = self.execution_client.start_execution(job=job)
        return execution_meta

    def process_executing_jobs(self):
        job_execution_states = self.get_job_execution_states()
        executed_jobs = [
            job for job in self.executing_jobs.values()
            if not self.job_is_executing(
                job=job, job_execution_states=job_execution_states)
        ]
        for executed_job in executed_jobs:
            self.process_executed_job(job=executed_job)

    def get_job_execution_states(self):
        job_execution_states = {}
        for job_key, job in self.executing_jobs.items():
            job_execution_state = self.execution_client.get_execution_state(
                job=job)
            job_execution_states[job_key] = job_execution_state
        return job_execution_states

    def job_is_executing(self, job=None, job_execution_states=None):
        if not job_execution_states.get(job['key']):
            is_executing = False
        else:
            is_executing = job_execution_states[job['key']]['executing']
        return is_executing

    def process_executed_job(self, job=None):
        raise NotImplementedError

    def complete_job(self, job=None):
        self.update_job(job=job['job'], updates={
            'status': self.job_client.Statuses.COMPLETED.name,
        })
        del self.job[job['key']]

    def update_job(self, job=None, updates=None):
        self.job_client.update_jobs(updates_by_uuid={
            job['uuid']: updates})

