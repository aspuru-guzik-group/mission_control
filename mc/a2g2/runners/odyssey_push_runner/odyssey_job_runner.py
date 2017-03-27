import collections
import types

from mc.job_runners.base_job_runner import BaseJobRunner
from mc.task_handlers.jobs.execute_job_task_handler import ExecuteJobTaskHandler
from mc.execution_clients.remote_slurm_execution_client import (
    RemoteSlurmExecutionClient as ExecutionClient)


class OdysseyJobRunner(object):
    def __init__(self, job_submission_factory=None,
                 ssh_client=None, **job_runner_kwargs):
        self.job_submission_factory = job_submission_factory
        self.ssh_client = ssh_client
        self.tasks_cfg = self.generate_tasks_cfg()
        self.base_job_runner = BaseJobRunner(
            task_handler=self.generate_task_handler(),
            default_job_tasks=self.tasks_cfg['default_job_tasks'],
            **job_runner_kwargs,
        )

    def generate_tasks_cfg(self):
        task_defs = collections.OrderedDict()
        task_defs['build_job'] = {
            'task_key': 'build_job',
            'task_type': 'build_job',
            'tick_fn': self.tick_build_job_task,
        }
        task_defs['execute_job'] = {
            'task_key': 'execute_job',
            'task_type': 'execute_job',
            'tick_fn': self.tick_execute_job_task,
        }
        task_defs['expose_outputs'] = {
            'task_key': 'expose_outputs',
            'task_type': 'expose_outputs',
            'tick_fn': self.tick_expose_outputs_task,
        }
        tasks_cfg = {
            'task_defs': task_defs,
            'default_job_tasks': [
                {k: v for k, v in task_def.items() if k != 'tick_fn'}
                for task_def in task_defs.values()
            ]
        }
        return tasks_cfg

    def generate_task_handler(self):
        task_handler = types.SimpleNamespace()
        task_handler.tick_task = self.tick_task
        return task_handler

    def tick_task(self, *args, task=None, **kwargs):
        task_type = task['task_type']
        task_def = self.tasks_cfg['task_defs'][task_type]
        task_def['tick_fn'](*args, task=task, **kwargs)

    def tick_build_job_task(self, *args, task=None, task_context=None,
                            **kwargs):
        job = task_context['job']
        submission = self.job_submission_factory.build_job_submission(job=job)
        execute_job_task = task_context['tasks']['execute_job']
        execute_job_task.setdefault('task_params', {})
        execute_job_task['task_params']['submission'] = submission
        task['status'] = 'COMPLETED'

    def tick_execute_job_task(self, *args, **kwargs):
        execution_client = ExecutionClient(ssh_client=self.ssh_client)
        task_handler = ExecuteJobTaskHandler(execution_client=execution_client)
        task_handler.tick_task(*args, **kwargs)

    def tick_expose_outputs_task(self, *args, task=None, task_context=None,
                                 **kwargs):
        job = task_context['job']
        job['data']['outputs'] = \
                task_context['tasks']['submit_job']['data']['outputs']

    def tick(self, *args, **kwargs):
        return self.base_job_runner.tick(*args, **kwargs)
