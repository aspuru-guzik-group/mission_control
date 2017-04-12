from mc.mc_utils.context_value_setter import set_context_values
from .base_task_handler import BaseTaskHandler


class SetValuesTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        set_context_values(value_specs=task['task_params']['value_specs'],
                           context=task_context)
        task['status'] = 'COMPLETED'
