import collections
import copy
import logging
import re

from mc.utils import dot_spec_loader

from .base_task_handler import BaseTaskHandler
from .wire_task_handler import WireTaskHandler
from .log_task_handler import LogTaskHandler


class McDefaultTaskHandler(object):
    def __init__(self, task_type_to_handler_dot_spec_fn=None, logger=None):
        self.logger = logger or logging
        self.task_type_to_handler_dot_spec = task_type_to_handler_dot_spec_fn \
                or self.default_task_type_to_handler_dot_spec

    def default_task_type_to_handler_dot_spec(self, task_type=None):
        handler_dot_spec = task_type
        if ':' not in handler_dot_spec:
            handler_dot_spec = re.sub(r'task(?!_handler)', r'task_handler',
                                      task_type)
            if not handler_dot_spec.endswith('task_handler'):
                handler_dot_spec += '_task_handler'
            handler_dot_spec += ':TaskHandler'
        return handler_dot_spec

    @classmethod
    def tick_task(cls, *args, **kwargs):
        return cls()._tick_task(*args, **kwargs)

    def _tick_task(self, *args, task_ctx=None, **kwargs):
        interpolated_task_ctx = self.get_interpolated_task_ctx(
            task_ctx=task_ctx)
        handler = self.get_handler_for_task_ctx(task_ctx=interpolated_task_ctx)
        handler.tick_task(*args, **kwargs)
        self.propagate_task_ctx_changes(
            orig_task_ctx=task_ctx, interpolated_task_ctx=interpolated_task_ctx)

    def get_handler_for_task_ctx(self, task_ctx=None):
        handler_cls = self.get_handler_cls_for_task_ctx(task_ctx=task_ctx)
        return handler_cls(task_ctx=task_ctx)

    def get_handler_cls_for_task_ctx(self, task_ctx=None):
        task_type = task_ctx['task']['task_type']
        if task_type == 'log': handler_cls = LogTaskHandler
        elif task_type == 'noop': handler_cls = NoOpTaskHandler
        elif task_type == 'print': handler_cls = PrintTaskHandler
        elif task_type == 'wire': handler_cls = WireTaskHandler
        else:
            handler_dot_spec = self.task_type_to_handler_dot_spec(task_type)
            handler_cls = self.load_from_dot_spec(dot_spec=handler_dot_spec)
        return handler_cls

    def get_interpolated_task_ctx(self, task_ctx=None):
        interpolated_task_ctx = {
            **task_ctx,
            'task': self.get_interpolated_task(task=task_ctx['task'],
                                               task_ctx=task_ctx)
        }
        return interpolated_task_ctx

    def get_interpolated_task(self, task=None, task_ctx=None):
        interpolated_task = copy.deepcopy(task)
        for parent, key, value in self.traverse_obj(interpolated_task):
            if self.value_has_interpolations(value=value):
                parent[key] = self.interpolate_value(value=value, ctx=task_ctx)
        return interpolated_task

    def traverse_obj(self, obj=None):
        parent = obj
        for key, value in self.get_kv_iterator_for_obj(obj=obj):
            if self.is_traversable(obj=value):
                yield from self.traverse_obj(obj=value)
            else: yield parent, key, value

    def is_traversable(self, obj=None):
        if isinstance(obj, str): return False
        return (isinstance(obj, collections.abc.Mapping)
                or isinstance(obj, collections.abc.Sequence))

    def get_kv_iterator_for_obj(self, obj=None):
        if isinstance(obj, collections.abc.Mapping): return obj.items()
        elif isinstance(obj, collections.abc.Sequence): return enumerate(obj)
        else: raise Exception("Could not get iterator for obj")

    def value_has_interpolations(self, value=None):
        return (
            isinstance(value, str)
            and (value.startswith('$') or value.startswith('\$'))
        )

    def interpolate_value(self, value=None, ctx=None):
        if value.startswith('$'):
            interpolated = self.get_ctx_value(ctx=ctx, dot_spec=value[1:])
        elif value.startswith('\$'): interpolated = value[1:]
        else: raise Exception("invalid value for interpolation")
        return interpolated

    def get_ctx_value(self, ctx=None, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.get_obj_value_from_dot_spec(
            obj={'ctx': ctx}, dot_spec=dot_spec)

    def load_from_dot_spec(self, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.load_from_dot_spec(
            dot_spec=dot_spec)

    def propagate_task_ctx_changes(self, orig_task_ctx=None,
                                   interpolated_task_ctx=None):
        for k, v in interpolated_task_ctx['task'].items():
            if k != 'task_params': orig_task_ctx['task'][k] = v

class NoOpTaskHandler(BaseTaskHandler):
    def tick_task(self, **kwargs): self.task['status'] = 'COMPLETED'

class PrintTaskHandler(BaseTaskHandler):
    def tick_task(self, **kwargs):
        print(self.task['task_params']['msg'])
        self.task['status'] = 'COMPLETED'
