import logging
import time


class BaseJobRunner(object):
    def __init__(self, job_spec_client=None, job_dir_factory=None,
                 execution_client=None, tick_interval=120, max_executing_jobs=3,
                 transfer_client=None, logger=None):
        self.job_spec_client = job_spec_client
        self.job_dir_factory = job_dir_factory
        self.execution_client = execution_client
        self.tick_interval = tick_interval
        self.max_executing_jobs = max_executing_jobs
        self.transfer_client = transfer_client
        self.logger = logger or logging

        self._ticking = False
        self.tick_counter = 0
        self.executing_jobs = {}
        self.transferring_jobs = {}

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
        self.process_transferring_jobs()
        self.process_executing_jobs()
        num_job_slots = self.max_executing_jobs - len(self.executing_jobs)
        if num_job_slots <= 0: return
        claimable_job_specs = self.fetch_claimable_job_specs()
        for claimable_job_spec in claimable_job_specs[:num_job_slots]:
            self.process_claimable_job_spec(job_spec=claimable_job_spec)

    def fetch_claimable_job_specs(self):
        logging.debug('fetch_claimable_job_specs')
        return self.job_spec_client.fetch_claimable_job_specs()

    def process_claimable_job_spec(self, job_spec=None):
        logging.debug('process_claimable_job_spec')
        claimed_spec = self.claim_job_spec(job_spec)
        if not claimed_spec: return
        dir_meta = self.build_job_dir(job_spec=claimed_spec)
        partial_job = {'key': claimed_spec['uuid'],
                       'job_spec': claimed_spec,
                       'dir': dir_meta}
        try:
            execution_meta = self.start_job_execution(job=partial_job)
            full_job = {**partial_job, 'execution': execution_meta}
            self.executing_jobs[partial_job['key']] = full_job
        except Exception as exception:
            logging.exception(exception)
            self.update_job_spec(job_spec=claimed_spec,
                                 updates={'status': 'Failed'})

    def claim_job_spec(self, job_spec=None):
        logging.debug('claim_job_spec')
        claimed_specs = self.job_spec_client.claim_job_specs(
            uuids=[job_spec['uuid']])
        return claimed_specs.get(job_spec['uuid'], False)

    def build_job_dir(self, job_spec=None):
        logging.debug('build_job_dir')
        job_dir_meta = self.job_dir_factory.build_dir_for_spec(
            job_spec=job_spec)
        return job_dir_meta

    def start_job_execution(self, job=None):
        logging.debug('start_job_execution')
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
        transfer_meta = self.start_job_transfer(job=job)
        self.transferring_jobs[job['key']] = {**job, 'transfer': transfer_meta}
        del self.executing_jobs[job['key']]

    def start_job_transfer(self, job=None):
        transfer_meta = self.transfer_client.start_transfer(job=job)
        return transfer_meta

    def process_transferring_jobs(self):
        job_transfer_states = self.get_job_transfer_states()
        transferred_jobs = [
            job for job in self.transferring_jobs.values()
            if not self.job_is_transferring(
                job=job, job_transfer_states=job_transfer_states)
        ]
        for transferred_job in transferred_jobs:
            self.process_transferred_job(job=transferred_job)

    def get_job_transfer_states(self):
        job_transfer_states = {}
        for job_key, job in self.transferring_jobs.items():
            job_transfer_state = self.transfer_client.get_transfer_state(
                job=job)
            job_transfer_states[job_key] = job_transfer_state
        return job_transfer_states

    def job_is_transferring(self, job=None, job_transfer_states=None):
        if not job_transfer_states.get(job['key']):
            is_transferring = False
        else:
            is_transferring = job_transfer_states[job['key']]['transferring']
        return is_transferring

    def process_transferred_job(self, job=None):
        logging.debug('process_transferred_job')
        self.update_job_spec(job_spec=job['job_spec'], updates={
            'status': self.job_spec_client.Statuses.Completed.name,
            'transfer_meta': job['transfer'],
        })
        del self.transferring_jobs[job['key']]

    def update_job_spec(self, job_spec=None, updates=None):
        self.job_spec_client.update_job_specs(updates_by_uuid={
            job_spec['uuid']: updates})

