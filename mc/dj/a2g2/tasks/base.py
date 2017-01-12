from uuid import uuid4
import logging

class BaseTask(object):
    def __init__(self, *args, uuid=None, key=None, status='PENDING', data=None,
                 logger=None, flow=None, input=None, output=None, error=None,
                 **kwargs):
        if uuid is None: uuid = str(uuid4())
        self.key = key
        self.status = status
        self.data = data or {}
        self.logger = logger or logging
        self.flow = flow
        self.input = input
        self.output = output
        self.error = error

    def tick(self, *args, **kwargs): raise NotImplementedError

    def increment_tick_counter(self):
        if 'ticks' not in self.data: self.data['ticks'] = 0
        self.data['ticks'] += 1

    def mark_as_failed(self, error=None):
        self.status = 'FAILED'
        self.error = error