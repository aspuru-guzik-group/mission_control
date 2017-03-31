import logging
from uuid import uuid4


class BaseTaskRunner(object):
    def __init__(self, get_tasks=None, get_task_context=None, task_handler=None,
                 logger=None):
        self.get_tasks = get_tasks or (lambda: [])
        self.get_task_context = get_task_context or (lambda: {})
        self.task_handler = task_handler
        self.logger = logger or logging

    def tick_tasks(self):
        try:
            tasks_status = None
            current_task = self.get_first_incomplete_task()
            if not current_task: tasks_status = 'COMPLETED'
            while current_task:
                self.tick_task(task=current_task)
                next_task = self.get_first_incomplete_task()
                if next_task is None: # no more tasks.
                    tasks_status = 'COMPLETED'
                    break
                elif next_task is current_task: # task didn't finish yet.
                    tasks_status = 'RUNNING'
                    break
                current_task = next_task
            return tasks_status
        except Exception as error:
            raise Exception("Error ticking tasks: {}".format(error))

    def get_first_incomplete_task(self):
        for task in self.get_tasks():
            if task.get('status', None) != 'COMPLETED': return task
        return None

    def tick_task(self, task=None):
        try:
            self.task_handler.tick_task(
                task=task,
                task_context=self.get_decorated_task_context()
            )
            if task.get('status', None) == 'FAILED':
                error = "Task with key '{task_key}' failed: {task_error}".format(
                    task_key=task.get('task_key', '<unknown>'),
                    task_error=task['data'].get('error', 'unknown error'))
                raise Exception(error)
        except Exception as exception:
            msg = "Failed to tick task with key '{key}'".format(
                key=task.get('task_key', '<unknown>'))
            self.logger.debug('task: {}'.format(task))
            self.logger.exception(msg)
            raise exception

    def get_decorated_task_context(self):
        decorated_task_context = {
            'tasks': self.get_keyed_tasks(tasks=self.get_tasks()),
            **self.get_task_context()
        }
        return decorated_task_context

    def get_keyed_tasks(self, tasks=None):
        keyed_tasks = {}
        for task in tasks:
            key = task.get('task_key', str(uuid4()))
            keyed_tasks[key] = task
        return keyed_tasks
