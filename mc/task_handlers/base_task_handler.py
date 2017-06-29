import abc
import logging
import traceback

from . import constants

class BaseTaskHandler(object):
    """An abstract base class to help with common task_handler logic.

    Provides utility methods:
       - tick_task class method: instantiates class, sets self.task_ctx
       - validate_task_ctx
       - validate_task_params

    Provides shortcuts:
       - self.task => task_ctx['task']

    Subclasses should at least implement initial_tick method.
    
    """
    TaskStatuses = constants.TaskStatuses

    class InvalidTaskError(Exception):
        def __init__(self, msg=None, task=None, **kwargs):
            msg = (msg or '') + "\ntask: '{task}'".format(task=task)
            super().__init__(msg, **kwargs)

    class InvalidTaskCtxError(InvalidTaskError): pass
    class InvalidTaskParamsError(InvalidTaskError): pass

    def __init__(self, task_ctx=None, logger=None, **kwargs):
        self.task_ctx = task_ctx
        self.logger = logger or logging

    @property
    def task(self): return self.task_ctx['task']

    @task.setter
    def task(self, value): self.task_ctx['task'] = value

    @classmethod
    def tick_task(cls, task_ctx=None, logger=None, **kwargs):
        cls(task_ctx=task_ctx, logger=logger, **kwargs)._tick_task(**kwargs)

    def _tick_task(self, **kwargs):
        try:
            self._ensure_task()
            try: self.validate_task_ctx()
            except Exception as exc: raise self.InvalidTaskCtxError() from exc
            self.increment_tick_counter()
            if self.task['data']['_tick_counter'] == 1:
                self.task['status'] = 'RUNNING'
                try: self.validate_task_params()
                except Exception as exc:
                    raise self.InvalidTaskParamsError() from exc
                self.initial_tick(**kwargs) 
            else: self.intermediate_tick(**kwargs)
        except Exception as exception:
            self.logger.exception(self.__class__.__name__ +".tick_task")
            self.task['status'] = 'FAILED'
            self.task['data']['error'] = traceback.format_exc()

    def _ensure_task(self):
        self.task.setdefault('data', {})
        self.task['data'].setdefault('_tick_counter', 0)

    def validate_task_ctx(self): pass

    def increment_tick_counter(self): self.task['data']['_tick_counter'] += 1

    def validate_task_params(self): pass

    @abc.abstractmethod
    def initial_tick(self, task=None, **kwargs): raise NotImplementedError

    def intermediate_tick(self, task=None, **kwargs): raise NotImplementedError
