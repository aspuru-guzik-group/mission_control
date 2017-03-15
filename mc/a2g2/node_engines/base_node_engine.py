import logging

class BaseNodeEngine(object):
    def __init__(self, logger=None, **kwargs):
        self.logger = logger or logging

    def tick_node(self, *args, node=None, ctx=None, **kwargs):
        raise NotImplementedError

    def increment_node_tick_counter(self, node=None):
        node.setdefault('data', {})
        node['data'].setdefault('ticks', 0)
        node['data']['ticks'] += 1

    def mark_node_as_failed(self, node=None, error=None):
        node['status'] = 'FAILED'
        node['error'] = error
