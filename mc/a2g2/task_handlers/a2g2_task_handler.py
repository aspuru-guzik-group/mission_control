from mc.mc_utils import dot_spec_loader

class A2G2TaskHandler(object):
    @classmethod
    def tick_task(cls, *args, task=None, **kwargs):
        task_type = task['task_type']
        if task_type == 'set_value':
            handler_dot_spec = 'mc.task_handlers.set_value_task_handler'
        else:
            handler_dot_spec = task_type.replace('task', 'task_handler')
        if ':' not in handler_dot_spec: handler_dot_spec += ':TaskHandler'
        handler_cls = dot_spec_loader.DotSpecLoader.load_from_dot_spec(
            dot_spec=handler_dot_spec)
        return handler_cls().tick_task(*args, task=task, **kwargs)
