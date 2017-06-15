import logging

from mc.utils import dot_spec_loader

from .base_task_handler import BaseTaskHandler
from .wire_task_handler import WireTaskHandler
from .log_task_handler import LogTaskHandler


class ModuleDispatchTaskHandler(object):

    def __init__(self, task_type_to_handler_dot_spec_fn=None, logger=None):
        self.logger = logger or logging
        self.task_type_to_handler_dot_spec = task_type_to_handler_dot_spec_fn \
                or self.default_task_type_to_handler_dot_spec

    def default_task_type_to_handler_dot_spec(self, task_type=None):
        handler_dot_spec = task_type
        if ':' not in handler_dot_spec:
            handler_dot_spec = task_type.replace('task', 'task_handler')
            if not handler_dot_spec.endswith('_task_handler'):
                handler_dot_spec += '_task_handler'
            handler_dot_spec += ':TaskHandler'
        return handler_dot_spec

    @classmethod
    def tick_task(cls, *args, **kwargs):
        return cls()._tick_task(*args, **kwargs)

    def _tick_task(self, *args, task_ctx=None, **kwargs):
        interpolated_task_ctx = self.get_interpolated_task_ctx(
            task_ctx=task_ctx)
        handler = self.get_handler(task_ctx=interpolated_task_ctx)
        handler.tick_task(*args, **kwargs)
        self.propagate_task_ctx_changes(
            orig_task_ctx=task_ctx, interpolated_task_ctx=interpolated_task_ctx)

    def get_handler(self, task_ctx=None):
        task_type = task_ctx['task']['task_type']
        if task_type == 'log': handler_cls = LogTaskHandler
        elif task_type == 'noop': handler_cls = NoOpTaskHandler
        elif task_type == 'print': handler_cls = PrintTaskHandler
        elif task_type == 'wire': handler_cls = WireTaskHandler
        else:
            handler_dot_spec = self.task_type_to_handler_dot_spec(task_type)
            if ':' not in handler_dot_spec: handler_dot_spec += ':TaskHandler'
            handler_cls = dot_spec_loader.DotSpecLoader.load_from_dot_spec(
                dot_spec=handler_dot_spec)
        return handler_cls(task_ctx=task_ctx)

    def get_interpolated_task_ctx(self, task_ctx=None):
        interpolated_task_ctx = {
            **task_ctx,
            'task': self.get_interpolated_task(task=task_ctx['task'],
                                               task_ctx=task_ctx)
        }
        return interpolated_task_ctx

    def get_interpolated_task(self, task=None, task_ctx=None):
        interpolated_task = {
            **task,
            'task_params': self.interpolate_task_params(
                task=task, task_ctx=task_ctx)
        }
        return interpolated_task

    def interpolate_task_params(self, task=None, task_ctx=None):
        interpolated_task_params = {
            key: self.interpolate_task_param_value(value=value,
                                                   task=task,
                                                   task_ctx=task_ctx)
            for key, value in task.get('task_params', {}).items()
        }
        return interpolated_task_params

    def interpolate_task_param_value(self, value=None, task=None,
                                     task_ctx=None):
        if isinstance(value, str) and value.startswith('_ctx:'):
            dot_spec = value.split('_ctx:')[1]
            interpolated = self.get_ctx_value(ctx=task_ctx,
                                              dot_spec=dot_spec)
        else: interpolated = value
        return interpolated

    def get_ctx_value(self, ctx=None, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.get_obj_value_from_dot_spec(
            obj=ctx, dot_spec=dot_spec)

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
