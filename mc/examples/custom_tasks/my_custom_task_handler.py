from mc.task_handlers.base_task_handler import BaseTaskHandler


class MyCustomTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        print("MyCustomTaskHandler.initial_tick")

    def intermediate_tick(self, task=None, task_context=None, **kwargs):
        print("MyCustomTaskHandler.intermediate_tick")
        task['status'] = 'COMPLETED'

TaskHandler = MyCustomTaskHandler
