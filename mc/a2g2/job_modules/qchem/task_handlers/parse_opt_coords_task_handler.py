from mc.task_handlers.base_task_handler import BaseTaskHandler


class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        task['data']['atoms'] = [1, 2, 3]
        task['status'] = 'COMPLETED'
