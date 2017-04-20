from mc.mc_utils.context_value_setter import set_context_value
from .base_task_handler import BaseTaskHandler


class SetValueTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        set_context_value(value_spec=task['task_params'], context=task_context)
        task['status'] = 'COMPLETED'

TaskHandler = SetValueTaskHandler
