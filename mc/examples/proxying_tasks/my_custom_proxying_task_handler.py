from mc.task_handlers.base_task_handler import BaseTaskHandler


class MyCustomProxyingTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        print("MyCustomTaskHandler.initial_tick")
        task['proxied_task'] = {
            'task_type': 'print',
            'task_params': {'msg': 'I am the proxied task'}
        }

TaskHandler = MyCustomProxyingTaskHandler
