import logging

class BaseNode(object):
    def __init__(self, *args, status=None, data=None, logger=None, **kwargs):
        self.status = status
        self.data = data or {}
        self.logger = logger or logging

    def tick(self): raise NotImplementedError
