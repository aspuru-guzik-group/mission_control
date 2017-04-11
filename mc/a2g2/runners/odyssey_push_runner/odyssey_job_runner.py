import collections
import copy
import json
import os
import tempfile
import types

from mc.job_runners.base_job_runner import BaseJobRunner
from mc.task_handlers.jobs.execute_job_task_handler import ExecuteJobTaskHandler
from .odyssey_execution_client import OdysseyExecutionClient


class OdysseyJobRunner(object):
    def __init__(self, job_submission_factory=None,
                 ssh_client=None, **job_runner_kwargs):
        self.job_submission_factory = job_submission_factory
        self.execution_client = OdysseyExecutionClient(ssh_client=ssh_client)
        self.tasks_cfg = self.generate_tasks_cfg()
        self.base_job_runner = BaseJobRunner(
            task_handler=self.generate_task_handler(),
            get_default_job_tasks=self.get_default_job_tasks,
            **job_runner_kwargs,
        )

    def generate_tasks_cfg(self):
        task_defs = collections.OrderedDict()
        task_defs['build_job_submission'] = {
            'task_key': 'build_job_submission',
            'task_type': 'build_job_submission',
            'tick_fn': self.tick_build_job_submission_task,
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

    def get_default_job_tasks(self, *args, **kwargs):
        return copy.deepcopy(self.tasks_cfg['default_job_tasks'])

    def generate_task_handler(self):
        task_handler = types.SimpleNamespace()
        task_handler.tick_task = self.tick_task
        return task_handler

    def tick_task(self, *args, task=None, **kwargs):
        task_type = task['task_type']
        task_def = self.tasks_cfg['task_defs'][task_type]
        task_def['tick_fn'](*args, task=task, **kwargs)

    def tick_build_job_submission_task(self, *args, task=None,
                                       task_context=None, **kwargs):
        job = task_context['job']
        submission_dir = tempfile.mkdtemp(prefix='sf.')
        self.prepare_job_inputs(job=job, submission_dir=submission_dir)
        submission = self.job_submission_factory.build_job_submission(
            job=job, submission_dir=submission_dir)
        execute_job_task = task_context['tasks']['execute_job']
        execute_job_task.setdefault('task_params', {})
        execute_job_task['task_params']['submission'] = submission
        task['status'] = 'COMPLETED'

    def prepare_job_inputs(self, job=None, submission_dir=None):
        inputs_dir = os.path.join(submission_dir, 'inputs')
        os.makedirs(inputs_dir, exist_ok=True)
        serialized_artifacts = job['job_spec'].get('inputs', {}).get(
            'serialized_artifacts', {})
        for artifact_key, serialized_artifact in serialized_artifacts.items():
            self.prepare_input_artifact(
                artifact_key=artifact_key,
                artifact=json.loads(serialized_artifact),
                inputs_dir=inputs_dir)

    def prepare_input_artifact(self, artifact_key=None, artifact=None,
                               inputs_dir=None):
        if artifact['artifact_type'] == 'a2g2.artifacts.odyssey':
            os.symlink(artifact['artifact_params']['path'],
                       os.path.join(inputs_dir, artifact_key))

    def tick_execute_job_task(self, *args, **kwargs):
        task_handler = ExecuteJobTaskHandler(
            execution_client=self.execution_client)
        task_handler.tick_task(*args, **kwargs)

    def tick_expose_outputs_task(self, *args, task=None, task_context=None,
                                 **kwargs):
        job = task_context['job']
        job['data']['artifact'] = \
                task_context['tasks']['execute_job']['data']['artifact']
        task['status'] = 'COMPLETED'

    def tick(self, *args, **kwargs):
        return self.base_job_runner.tick(*args, **kwargs)
