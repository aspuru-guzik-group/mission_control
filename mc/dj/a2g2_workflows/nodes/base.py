from uuid import uuid4
import logging

class BaseNode(object):
    def __init__(self, *args, id=None, status=None, data=None, logger=None,
                 **kwargs):
        if id is None: id = uuid4()
        self.id = id
        self.status = status
        self.data = data or {}
        self.logger = logger or logging

    def tick(self, *args, **kwargs): raise NotImplementedError
