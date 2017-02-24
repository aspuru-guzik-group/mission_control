import logging
import tempfile
import time


class BaseJobRunner(object):
    def __init__(self, job_client=None, job_submission_factory=None,
                 execution_client=None, action_processor=None,
                 tick_interval=120, max_executing_jobs=3, 
                 get_job_execution_result=None, logger=None):
        self.job_client = job_client
        self.job_submission_factory = job_submission_factory
        self.execution_client = execution_client
        self.action_processor = action_processor
        self.tick_interval = tick_interval
        self.max_executing_jobs = max_executing_jobs
        self._get_job_execution_result = get_job_execution_result
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
        self.logger.debug('%s, tick #%s' % (self, self.tick_counter))
        self.process_executing_jobs()
        num_job_slots = self.max_executing_jobs - len(self.executing_jobs)
        if num_job_slots <= 0: return
        claimable_jobs = self.fetch_claimable_jobs()
        for claimable_job in claimable_jobs[:num_job_slots]:
            self.process_claimable_job(job=claimable_job)

    def fetch_claimable_jobs(self):
        self.logger.debug('fetch_claimable_jobs')
        return self.job_client.fetch_claimable_jobs()

    def process_claimable_job(self, job=None):
        self.logger.debug('process_claimable_job')
        claimed_job = self.claim_job(job)
        if not claimed_job: return
        self.register_job(job=claimed_job)
        try:
            try:
                claimed_job['submission'] = self.build_job_submission(
                    job=claimed_job)
            except Exception as error:
                raise Exception("Error building job directory: %s" % error)
            try:
                claimed_job['execution'] = self.start_job_execution(
                    job=claimed_job)
            except Exception as error:
                raise Exception("Error starting job execution: %s" % error)
        except Exception as error:
            self.fail_job(job=claimed_job, error=error)

    def register_job(self, job=None):
        self.jobs[job['uuid']] = job

    def unregister_job(self, job=None):
        for registry in [self.jobs, self.executing_jobs]:
            if job['uuid'] in registry: del registry[job['uuid']]

    def fail_job(self, job=None, error=None):
        self.logger.exception(error)
        self.update_job(job=job, updates={'status': 'FAILED',
                                          'error': str(error)})
        self.unregister_job(job=job)

    def claim_job(self, job=None):
        self.logger.debug('claim_job')
        claimed_jobs = self.job_client.claim_jobs(
            uuids=[job['uuid']])
        return claimed_jobs.get(job['uuid'], False)

    def build_job_submission(self, job=None):
        self.logger.debug('build_job_submission')
        job['job_dir'] = self.generate_job_dir(job=job)
        actions =  job.get('job_spec', {}).get('pre_build_actions', None)
        if actions: self.process_actions(actions=actions, job=job)
        job_submission_meta = self.job_submission_factory.build_job_submission(
            job=job, output_dir=job['job_dir'])
        return job_submission_meta

    def generate_job_dir(self, job=None):
        return tempfile.mkdtemp(prefix='job_dir.{}'.format(job['uuid']))

    def process_actions(self, actions=None, job=None):
        for action in actions:
            self.action_processor.process_action(action=action, ctx=job)

    def start_job_execution(self, job=None):
        self.logger.debug('start_job_execution, %s' % job)
        self.executing_jobs[job['uuid']] = job
        execution_meta = self.execution_client.start_execution(job=job)
        return execution_meta

    def process_executing_jobs(self):
        self.logger.debug('process_executing_jobs')
        job_execution_states = self.get_job_execution_states()
        executed_jobs = [
            job for job in self.executing_jobs.values()
            if not self.job_is_executing(
                job=job, job_execution_states=job_execution_states)
        ]
        for executed_job in executed_jobs:
            del self.executing_jobs[executed_job['uuid']]
            self.process_executed_job(job=executed_job)

    def get_job_execution_states(self):
        job_execution_states = {}
        for job_key, job in self.executing_jobs.items():
            job_execution_state = self.execution_client.get_execution_state(
                job=job)
            job_execution_states[job_key] = job_execution_state
        return job_execution_states

    def job_is_executing(self, job=None, job_execution_states=None):
        if not job_execution_states.get(job['uuid']):
            is_executing = False
        else:
            is_executing = job_execution_states[job['uuid']]['executing']
        return is_executing

    def process_executed_job(self, job=None):
        self.logger.debug('process_executed_job')
        job['completed_dir'] = job['execution']['completed_dir']
        job['execution']['result'] = self.get_job_execution_result(job=job)
        actions = job.get('job_spec', {}).get('post_exec_actions', None)
        if actions: self.process_actions(actions=actions, job=job)
        self.complete_job(job=job)

    def get_job_execution_result(self, job=None):
        _get_job_execution_result = getattr(self, '_get_job_execution_result',
                                            None)
        if _get_job_execution_result:
            result = _get_job_execution_result(job_state=job)
        else:
            result = 'COMPLETED'
        return result

    def complete_job(self, job=None):
        self.logger.debug('complete_job')
        execution_result = job['execution']['result']
        if execution_result['result'] == 'COMPLETED':
            self.update_job(job=job, updates={
                'status': 'COMPLETED',
                'data': job.get('data', {})
            })
            self.unregister_job(job=job)
        else:
            self.fail_job(
                job=job,
                error=execution_result.get('error', 'error unknown')
            )

    def update_job(self, job=None, updates=None):
        self.job_client.update_jobs(updates_by_uuid={
            job['uuid']: updates})
