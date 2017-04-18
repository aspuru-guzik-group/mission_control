import logging


class BaseTaskHandler(object):
    def __init__(self, *args, logger=None, **kwargs):
        self.logger = logger or logging

    def tick_task(self, *args, task=None, **kwargs):
        try:
            self._ensure_task(task=task)
            self.increment_tick_counter(task=task)
            if task['data']['_tick_counter'] == 1:
                task['status'] = 'RUNNING'
                self.initial_tick(task=task, **kwargs)
            else:
                self.intermediate_tick(task=task, **kwargs)
        except Exception as exception:
            self.logger.exception(self.__class__.__name__ +".tick_task")
            task['status'] = 'FAILED'
            task['data']['error'] = self.stringify_exception(exception)

    def stringify_exception(self, exception=None):
        return '[%s] %s)' % (type(exception), exception)

    def _ensure_task(self, task=None):
        task.setdefault('data', {})
        task['data'].setdefault('_tick_counter', 0)

    def increment_tick_counter(self, task=None):
        task['data']['_tick_counter'] += 1

    def initial_tick(self, task=None, **kwargs):
        raise NotImplementedError

    def intermediate_tick(self, task=None, **kwargs):
        raise NotImplementedError
