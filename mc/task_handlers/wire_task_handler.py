from mc.utils import context_value_setter
from .base_task_handler import BaseTaskHandler


class WireTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        for wiring in task['task_params'].get('wirings', []):
            self.process_wiring(wiring=wiring, task_context=task_context)
        task['status'] = 'COMPLETED'

    def process_wiring(self, wiring=None, task_context=None):
        context_value_setter.set_context_value(
            value_spec=wiring,
            context=task_context
        )

TaskHandler = WireTaskHandler
