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
        self.get_default_job_tasks = get_default_job_tasks or \
            self._get_default_job_tasks 
        self.logger = logger or logging

        self._ticking = False
        self.tick_counter = 0
        self.jobs = {}
        self.running_jobs = {}

    def _get_default_job_tasks(self, *args, **kwargs): return []

    def start(self):
        self._ticking = True
        self.run()

    def stop(self): self._ticking = False

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
        num_job_slots = self.max_running_jobs - len(self.running_jobs)
        if num_job_slots <= 0: return
        for job in self.claim_jobs(params={'limit': num_job_slots}):
            self.start_job(job=job)
        self.process_running_jobs()

    def job_is_running(self, job=None):
        return job['status'] == 'RUNNING'

    def claim_jobs(self, *args, **kwargs):
        return self.job_client.claim_jobs(*args, **kwargs)

    def start_job(self, job=None):
        try:
            self.register_job(job=job)
            if 'tasks' not in job:
                job['tasks'] = self.get_default_job_tasks(job=job)
            self.running_jobs[job['uuid']] = job
        except Exception as exception:
            self.logger.exception("Failed to start job")
            self.fail_job(job=job, error=traceback.format_exc())

    def register_job(self, job=None):
        self.jobs[job['uuid']] = job

    def fail_job(self, job=None, error=None):
        job.update({'status': 'FAILED', 'error': error})
        self.unregister_job(job=job)

    def unregister_job(self, job=None):
        for registry in [self.jobs, self.running_jobs]:
            if job['uuid'] in registry: del registry[job['uuid']]

    def process_running_jobs(self):
        keyed_running_jobs = list(self.running_jobs.items())
        for key, job in keyed_running_jobs:
            if self.job_is_running(job=job):
                try: self.tick_job(job=job)
                except Exception as error: self.fail_job(job=job, error=error)
            else: self.unregister_job(job=job)
        self.update_jobs(jobs=[job for key, job in keyed_running_jobs])

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

    def update_jobs(self, jobs=None):
        self.job_client.update_jobs(updates_by_uuid={
            job['uuid']: {
                attr: job.get(attr)
                for attr in ['data', 'status', 'error']
            } for job in jobs
        })
