import shlex

from mc.mc_utils import dot_spec_loader


class A2G2TaskHandler(object):
    @classmethod
    def tick_task(cls, *args, **kwargs):
        return cls()._tick_task(*args, **kwargs)

    def _tick_task(self, *args, task=None, task_context=None, **kwargs):
        handler = self.get_handler(task=task)
        interpolated_task = self.get_interpolated_task(
            task=task, task_context=task_context)
        result = handler.tick_task(
            *args,
            task=interpolated_task,
            task_context={**task_context, 'task': interpolated_task},
            **kwargs)
        for k, v in interpolated_task.items():
            if k != 'task_params': task[k] = v
        return result

    def get_handler(self, task=None):
        task_type = task['task_type']
        if task_type == 'noop': return NoOpTaskHandler()
        if task_type == 'print': return PrintTaskHandler()
        if task_type == 'set_value':
            handler_dot_spec = 'mc.task_handlers.set_value_task_handler'
        else:
            handler_dot_spec = task_type.replace('task', 'task_handler')
        if ':' not in handler_dot_spec: handler_dot_spec += ':TaskHandler'
        handler_cls = dot_spec_loader.DotSpecLoader.load_from_dot_spec(
            dot_spec=handler_dot_spec)
        return handler_cls()


    def get_interpolated_task(self, task=None, task_context=None):
        interpolated_task = {
            **task,
            'task_params': self.interpolate_task_params(
                task=task, task_context=task_context)
        }
        return interpolated_task

    def interpolate_task_params(self, task=None, task_context=None):
        interpolated_task_params = {
            key: self.interpolate_task_param_value(value=value,
                                                   task=task,
                                                   task_context=task_context)
            for key, value in task.get('task_params', {}).items()
        }
        return interpolated_task_params

    def interpolate_task_param_value(self, value=None, task=None,
                                     task_context=None):
        if isinstance(value, str) and value.startswith('_ctx:'):
            dot_spec = value.split('_ctx:')[1]
            interpolated = self.get_ctx_value(ctx=task_context,
                                              dot_spec=dot_spec)
        else:
            interpolated = value
        return interpolated

    def get_ctx_value(self, ctx=None, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.get_obj_value_from_dot_spec(
            obj=ctx, dot_spec=dot_spec)

    @classmethod
    def compile_tasks(cls, *args, **kwargs):
        return cls()._compile_tasks(*args, **kwargs)

    def _compile_tasks(self, tasks=None):
        if not tasks: return tasks
        return [self.compile_task(task=task) for task in tasks]

    def compile_task(self, task=None):
        if isinstance(task, str):
            compiled_task = self.compile_str_task(str_task=task)
        else: compiled_task = task
        return compiled_task

    def compile_str_task(self, str_task=None):
        task_type, *param_key_value_strs = shlex.split(str_task)
        compiled_task = {
            'task_type': task_type,
            'task_params': self.kv_strs_to_dict(kv_strs=param_key_value_strs)
        }
        return compiled_task

    def kv_strs_to_dict(self, kv_strs=None):
        return dict([kv_str.split('=', 1) for kv_str in kv_strs])

class NoOpTaskHandler(object):
    def tick_task(self, task=None, **kwargs):
        task['status'] = 'COMPLETED'

class PrintTaskHandler(object):
    def tick_task(self, task=None, **kwargs):
        print(task['task_params']['message'])
        task['status'] = 'COMPLETED'
