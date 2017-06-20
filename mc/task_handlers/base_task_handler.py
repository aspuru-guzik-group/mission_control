import logging
import traceback

from . import constants

class BaseTaskHandler(object):
    TaskStatuses = constants.TaskStatuses

    class InvalidTaskParamsError(Exception): pass

    def __init__(self, *args, task_ctx=None, logger=None, **kwargs):
        self.task_ctx = task_ctx
        self.logger = logger or logging

    @property
    def task(self): return self.task_ctx['task']

    @task.setter
    def task(self, value): self.task_ctx['task'] = value

    @classmethod
    def tick_task(cls, *args, task_ctx=None, logger=None, **kwargs):
        cls(*args, task_ctx=task_ctx, logger=logger, **kwargs)._tick_task(
            *args, **kwargs)

    def _tick_task(self, *args, **kwargs):
        try:
            self._ensure_task()
            self.increment_tick_counter()
            if self.task['data']['_tick_counter'] == 1:
                self.task['status'] = 'RUNNING'
                self.initial_tick(**kwargs)
            else: self.intermediate_tick(**kwargs)
        except Exception as exception:
            self.logger.exception(self.__class__.__name__ +".tick_task")
            self.task['status'] = 'FAILED'
            self.task['data']['error'] = traceback.format_exc()

    def _ensure_task(self):
        self.task.setdefault('data', {})
        self.task['data'].setdefault('_tick_counter', 0)

    def increment_tick_counter(self): self.task['data']['_tick_counter'] += 1

    def initial_tick(self, task=None, **kwargs): raise NotImplementedError

    def intermediate_tick(self, task=None, **kwargs): raise NotImplementedError
