from mc.mc_utils import dot_spec_loader

class DotSpecTaskHandler(object):
    @classmethod
    def tick_task(cls, *args, task=None, **kwargs):
        task_type = task['task_type']
        # hack to set convention for default task handler
        if ':' not in task_type: task_type = task_type + ':TaskHandler'
        handler_cls = dot_spec_loader.DotSpecLoader.load_from_dot_spec(
            dot_spec=task_type)
        return handler_cls().tick_task(*args, task=task, **kwargs)
