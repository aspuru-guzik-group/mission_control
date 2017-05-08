import logging

from mc.task_runners.base_task_runner import BaseTaskRunner


class JobTicker(object):
    def __init__(self, task_handler=None, logger=None):
        self.task_handler = task_handler
        self.logger = logger or logging

    def tick_job(self, job=None):
        task_runner = self.get_task_runner(job=job)
        tasks_status = task_runner.tick_tasks()
        job['status'] = tasks_status

    def get_task_runner(self, job=None):
        return BaseTaskRunner(
            get_tasks=lambda: job.get('tasks', []),
            get_task_context={'job': job},
            task_handler=self.task_handler
        )

