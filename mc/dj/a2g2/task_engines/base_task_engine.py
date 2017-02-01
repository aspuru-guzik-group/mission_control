import logging

class BaseTaskEngine(object):
    def __init__(self, logger=None, **kwargs):
        self.logger = logger or logging

    def tick_task(self, *args, task=None, ctx=None, **kwargs):
        raise NotImplementedError

    def increment_task_tick_counter(self, task=None):
        task.setdefault('data', {})
        task['data'].setdefault('ticks', 0)
        task['data']['ticks'] += 1

    def mark_task_as_failed(self, task=None, error=None):
        task['status'] = 'FAILED'
        task['error'] = error
