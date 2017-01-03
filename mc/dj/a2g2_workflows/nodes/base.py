from uuid import uuid4
import logging

class BaseNode(object):
    def __init__(self, *args, id=None, status='PENDING', data=None, logger=None,
                 workflow=None, **kwargs):
        if id is None: id = uuid4()
        self.id = id
        self.status = status
        self.data = data or {}
        self.logger = logger or logging
        self.workflow = workflow

    def tick(self, *args, **kwargs): raise NotImplementedError

    def increment_tick_counter(self):
        if 'ticks' not in self.data: self.data['ticks'] = 0
        self.data['ticks'] += 1
