from mc.task_handlers.base_proxying_task_handler import BaseProxyingTaskHandler


class MyCustomProxyingTaskHandler(BaseProxyingTaskHandler):
    def generate_proxied_task(self, *args, **kwargs):
        print("MyCustomTaskHandler.generate_proxying_task")
        return {
            'task_type': 'print',
            'task_params': {'msg': 'I am the proxied task'}
        }

TaskHandler = MyCustomProxyingTaskHandler
