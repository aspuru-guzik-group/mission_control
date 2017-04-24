import logging
import time
import traceback

from mc.task_runners.base_task_runner import BaseTaskRunner


class BaseJobRunner(object):
    def __init__(self, job_client=None, tick_interval=120, max_running_jobs=3, 
                 task_handler=None, get_default_job_tasks=None, logger=None):
        self.job_client = job_client
        self.tick_interval = tick_interval
        self.max_running_jobs = max_running_jobs
        self.task_handler = task_handler
        if get_default_job_tasks:
            self.get_default_job_tasks = get_default_job_tasks
        self.logger = logger or logging

        self._ticking = False
        self.tick_counter = 0
        self.jobs = {}
        self.running_jobs = {}

    def get_default_job_tasks(self, job=None): return []

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
        self.process_running_jobs()
        num_job_slots = self.max_running_jobs - len(self.running_jobs)
        if num_job_slots <= 0: return
        claimable_jobs = self.fetch_claimable_jobs()
        for claimable_job in claimable_jobs[:num_job_slots]:
            self.process_claimable_job(job=claimable_job)

    def process_running_jobs(self):
        completed_jobs = {}
        for key, job in list(self.running_jobs.items()):
            if self.job_is_running(job=job):
                try: self.tick_job(job=job)
                except Exception as error: self.fail_job(job=job, error=error)
            else:
                completed_jobs[key] = job
        for key, job in completed_jobs.items(): self.complete_job(job=job)

    def job_is_running(self, job=None):
        return job['status'] == 'RUNNING'

    def fetch_claimable_jobs(self):
        return self.job_client.fetch_claimable_jobs()

    def process_claimable_job(self, job=None):
        claimed_job = self.claim_job(job)
        if not claimed_job: return
        job = claimed_job
        self.register_job(job=job)
        try: self.start_job(job=job)
        except Exception as exception:
            self.logger.exception("Job failed")
            self.fail_job(job=job, error=traceback.format_exc())

    def claim_job(self, job=None):
        claimed_jobs = self.job_client.claim_jobs(
            uuids=[job['uuid']])
        return claimed_jobs.get(job['uuid'], False)

    def register_job(self, job=None):
        self.jobs[job['uuid']] = job

    def unregister_job(self, job=None):
        for registry in [self.jobs, self.running_jobs]:
            if job['uuid'] in registry: del registry[job['uuid']]

    def fail_job(self, job=None, error=None):
        self.update_job(job=job, updates={'status': 'FAILED',
                                          'error': error})
        self.unregister_job(job=job)

    def start_job(self, job=None):
        if 'tasks' not in job:
            job['tasks'] = self.get_default_job_tasks(job=job)
        self.running_jobs[job['uuid']] = job
        self.tick_job(job=job)

    def tick_job(self, job=None):
        try:
            task_runner = self.get_task_runner(job=job)
            tasks_status = task_runner.tick_tasks()
            job['status'] = tasks_status
        except Exception as exception:
            self.fail_job(job=job, error=traceback.format_exc())

    def get_task_runner(self, job=None):
        return BaseTaskRunner(
            get_tasks=lambda: job.get('tasks', []),
            get_task_context=lambda : self.get_task_context(job=job),
            task_handler=self.task_handler
        )

    def get_task_context(self, job=None):
        task_context = {'job': job}
        return task_context

    def complete_job(self, job=None):
        if job['status'] == 'FAILED':
            self.fail_job(job=job, error=job.get('error', '<unknown error>'))
        else:
            self.update_job(job=job, updates={
                'status': 'COMPLETED',
                'data': job.get('data', {})
            })
            self.unregister_job(job=job)

    def update_job(self, job=None, updates=None):
        self.job_client.update_jobs(updates_by_uuid={job['uuid']: updates})
