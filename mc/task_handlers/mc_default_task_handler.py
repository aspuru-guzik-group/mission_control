import collections
import copy
import re

from mc.utils import dot_spec_loader

from .base_task_handler import BaseTaskHandler
from .flow_task_handler import FlowTaskHandler
from .job_task_handler import JobTaskHandler
from .log_task_handler import LogTaskHandler
from .wire_task_handler import WireTaskHandler
from .spread_task_handler import SpreadTaskHandler
from .switch_task_handler import SwitchTaskHandler
from . import constants


class McDefaultTaskHandler(BaseTaskHandler):
    """Defult MissionControl TaskHandler.

    Dispatches to tick functions by switching on
    task_ctx['task']['task_type'].

    Also performs '$' interpolation for doing simple wirings.
    """
    NO_INTERPOLATE_KEY = constants.NO_INTERPOLATE_KEY

    class InvalidTaskTypeError(BaseTaskHandler.InvalidTaskError):
        pass

    def __init__(self, task_type_to_tick_fn_dot_spec_fn=None, **kwargs):
        super().__init__(**kwargs)
        self.task_type_to_tick_fn_dot_spec = (
            task_type_to_tick_fn_dot_spec_fn
            or self.default_task_type_to_tick_fn_dot_spec
        )

    def default_task_type_to_tick_fn_dot_spec(self, task_type=None):
        """@TODO

        """
        tick_fn_dot_spec = task_type
        if ':' not in tick_fn_dot_spec:
            tick_fn_dot_spec = re.sub(r'task(?!_handler)', r'task_handler',
                                      task_type)
            if not tick_fn_dot_spec.endswith('task_handler'):
                tick_fn_dot_spec += '_task_handler'
            tick_fn_dot_spec += ':TaskHandler.tick_task'
        return tick_fn_dot_spec

    def _tick_task(self, *args, **kwargs):
        if self.should_interpolate_task_ctx():
            self.interpolate_and_tick_task_ctx(*args, **kwargs)
        else:
            self.tick_task_ctx(*args, task_ctx=self.task_ctx, **kwargs)

    def interpolate_and_tick_task_ctx(self, *args, **kwargs):
        interpolated_task_ctx = self.get_interpolated_task_ctx()
        self.tick_task_ctx(*args, task_ctx=interpolated_task_ctx, **kwargs)
        self.task.update(interpolated_task_ctx['task'])

    def should_interpolate_task_ctx(self):
        if self.NO_INTERPOLATE_KEY in self.task:
            return False
        if self.get_task_type() == 'mc.tasks.flow':
            return False
        return True

    def get_task_type(self, task_ctx=None):
        task_ctx = task_ctx or self.task_ctx
        try:
            return task_ctx['task']['task_type']
        except Exception as exc:
            task = task_ctx.get('task')
            raise self.InvalidTaskError(
                msg="invalid task_type", task=task) from exc

    def get_interpolated_task_ctx(self):
        interpolated_task_ctx = {**self.task_ctx,
                                 'task': self.get_interpolated_task()}
        return interpolated_task_ctx

    def get_interpolated_task(self):
        interpolated_task = copy.deepcopy(self.task)
        for parent, key, value in self.traverse_obj(interpolated_task):
            if self.value_has_interpolations(value=value):
                parent[key] = self.interpolate_value(
                    value=value, ctx=self.task_ctx)
        interpolated_task[self.NO_INTERPOLATE_KEY] = True
        return interpolated_task

    def traverse_obj(self, obj=None, exclude_keys=None):
        exclude_keys = exclude_keys or {}
        parent = obj
        for key, value in self.get_kv_iterator_for_obj(obj=obj):
            if key in exclude_keys:
                continue
            if self.is_traversable(obj=value):
                yield from self.traverse_obj(
                    obj=value, exclude_keys=exclude_keys)
            else:
                yield parent, key, value

    def is_traversable(self, obj=None):
        if isinstance(obj, str):
            return False
        return (
            (isinstance(obj, collections.abc.Mapping)
             or isinstance(obj, collections.abc.Sequence))
        )

    def get_kv_iterator_for_obj(self, obj=None):
        if isinstance(obj, collections.abc.Mapping):
            return obj.items()
        elif isinstance(obj, collections.abc.Sequence):
            return enumerate(obj)
        else:
            raise Exception("Could not get iterator for obj")

    def value_has_interpolations(self, value=None):
        return (
            isinstance(value, str)
            and (value.startswith('$') or value.startswith('\$'))
        )

    def interpolate_value(self, value=None, ctx=None):
        if value.startswith('$'):
            interpolated = self.get_ctx_value(ctx=ctx, dot_spec=value[1:])
        elif value.startswith('\$'):
            interpolated = value[1:]
        else:
            raise Exception("invalid value for interpolation")
        return interpolated

    def get_ctx_value(self, ctx=None, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.get_obj_value_from_dot_spec(
            obj={'ctx': ctx}, dot_spec=dot_spec)

    def load_from_dot_spec(self, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.load_from_dot_spec(
            dot_spec=dot_spec)

    def tick_task_ctx(self, *args, task_ctx=None, **kwargs):
        tick_fn = self.get_tick_fn_for_task_ctx(task_ctx=task_ctx, **kwargs)
        tick_fn(*args, task_ctx=task_ctx, **kwargs)

    def get_tick_fn_for_task_ctx(self, task_ctx=None,
                                 tick_fn_for_task_ctx_overrides=None,
                                 **kwargs):
        task_type = self.get_task_type(task_ctx=task_ctx)
        if task_type in (tick_fn_for_task_ctx_overrides or {}):
            tick_fn = tick_fn_for_task_ctx_overrides[task_type]
        elif task_type == 'flow':
            tick_fn = FlowTaskHandler.tick_task
        elif task_type == 'job':
            tick_fn = JobTaskHandler.tick_task
        elif task_type == 'log':
            tick_fn = LogTaskHandler.tick_task
        elif task_type == 'noop':
            tick_fn = NoOpTaskHandler.tick_task
        elif task_type == 'print':
            tick_fn = PrintTaskHandler.tick_task
        elif task_type == 'spread':
            tick_fn = SpreadTaskHandler.tick_task
        elif task_type == 'switch':
            tick_fn = SwitchTaskHandler.tick_task
        elif task_type == 'wire':
            tick_fn = WireTaskHandler.tick_task
        else:
            tick_fn_dot_spec = self.task_type_to_tick_fn_dot_spec(task_type)
            tick_fn = self.load_from_dot_spec(dot_spec=tick_fn_dot_spec)
        return tick_fn


class NoOpTaskHandler(BaseTaskHandler):
    def initial_tick(self, *args, **kwargs):
        self.task['status'] = 'COMPLETED'


class PrintTaskHandler(BaseTaskHandler):
    def initial_tick(self, *args, **kwargs):
        print(self.task['task_params']['msg'])
        self.task['status'] = 'COMPLETED'
