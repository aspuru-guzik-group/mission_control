import logging
import tempfile
import time


class BaseJobRunner(object):
    def __init__(self, job_client=None, job_submission_factory=None,
                 execution_client=None, task_runner=None,
                 tick_interval=120, max_running_jobs=3, logger=None):
        self.job_client = job_client
        self.job_submission_factory = job_submission_factory
        self.execution_client = execution_client
        self.task_runner = task_runner
        self.tick_interval = tick_interval
        self.max_running_jobs = max_running_jobs
        self.logger = logger or logging

        self.default_tasks = [{'type': 'job:execute'}]

        self._ticking = False
        self.tick_counter = 0
        self.jobs = {}
        self.running_jobs = {}

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
                try: self.tick_job_tasks(job=job)
                except Exception as error: self.fail_job(job=job, error=error)
            else:
                completed_jobs[key] = job
        for key, job in completed_jobs.items(): self.complete_job(job=job)

    def job_is_running(self, job=None):
        incomplete_task = self.get_first_incomplete_task_for_job(job=job)
        if incomplete_task: return True
        return False

    def fetch_claimable_jobs(self):
        return self.job_client.fetch_claimable_jobs()

    def process_claimable_job(self, job=None):
        claimed_job = self.claim_job(job)
        if not claimed_job: return
        job = claimed_job
        self.register_job(job=job)
        try:
            job['submission'] = self.build_job_submission(job=job)
            self.start_job(job=job)
        except Exception as error:
            self.fail_job(job=job, error=error)

    def register_job(self, job=None):
        self.jobs[job['uuid']] = job

    def unregister_job(self, job=None):
        for registry in [self.jobs, self.running_jobs]:
            if job['uuid'] in registry: del registry[job['uuid']]

    def fail_job(self, job=None, error=None):
        self.logger.exception(error)
        self.update_job(job=job, updates={'status': 'FAILED',
                                          'error': str(error)})
        self.unregister_job(job=job)

    def claim_job(self, job=None):
        claimed_jobs = self.job_client.claim_jobs(
            uuids=[job['uuid']])
        return claimed_jobs.get(job['uuid'], False)

    def build_job_submission(self, job=None):
        try:
            job['job_dir'] = self.generate_job_dir(job=job)
            job_submission_meta = self.job_submission_factory\
                    .build_job_submission(job=job, output_dir=job['job_dir'])
            return job_submission_meta
        except Exception as error:
            raise Exception("Error building job submission: %s" % error)

    def generate_job_dir(self, job=None):
        return tempfile.mkdtemp(prefix='job_dir.{}'.format(job['uuid']))

    def start_job(self, job=None):
        if 'tasks' not in job: job['tasks'] = self.default_tasks
        self.running_jobs[job['uuid']] = job
        self.tick_job_tasks(job=job)

    def tick_job_tasks(self, job=None):
        try:
            current_task = self.get_first_incomplete_task_for_job(job=job)
            while current_task:
                self.tick_task(task=current_task, job=job)
                if current_task.get('status', None) == 'COMPLETED':
                    next_task = self.get_first_incomplete_task_for_job(job=job)
                    if next_task is current_task: next_task = None
                else: next_task = None
                current_task = next_task
        except Exception as error:
            raise Exception("Error ticking job tasks: {}".format(error))

    def get_first_incomplete_task_for_job(self, job=None):
        for task in job.get('tasks', []):
            if task.get('status', None) != 'COMPLETED': return task
        return None

    def tick_task(self, task=None, job=None):
        try:
            self.task_runner.tick_task(task=task, job=job)
            if task.get('status', None) == 'FAILED':
                error = "Task '{task}' failed: {task_error}".format(
                    task=task,
                    task_error=task['state'].get('error', 'unknown error'))
                raise Exception(error)
        except Exception as error:
            error = "Failed to tick task '{task}': {error}".format(
                task=task, error=error)
            raise Exception(error)

    def complete_job(self, job=None):
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
