from mc.task_handlers.base_task_handler import BaseTaskHandler


class MyCustomProxyingTaskHandler(BaseTaskHandler):
    def initial_tick(self, *args, **kwargs):
        print("MyCustomTaskHandler.initial_tick")
        self.task['proxied_task'] = {
            'task_type': 'print',
            'task_params': {'msg': 'I am the proxied task'}
        }

TaskHandler = MyCustomProxyingTaskHandler
