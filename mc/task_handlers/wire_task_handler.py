from mc.utils import context_value_setter
from .base_task_handler import BaseTaskHandler


class WireTaskHandler(BaseTaskHandler):
    """Performs wirings.

    Dispatches to :meth:`mc.utils.context_value_setter.set_context_value`

    task_params:
        wirings: a list of value_specs,
        per :meth:`mc.utils.context_value_setter.set_context_value`
    """

    def validate_task_params(self):
        assert self.task['task_params']['wirings']

    def initial_tick(self):
        for wiring in self.task['task_params'].get('wirings', []):
            self.process_wiring(wiring=wiring)
        self.task['status'] = 'COMPLETED'

    def process_wiring(self, wiring=None):
        context_value_setter.set_context_value(value_spec=wiring,
                                               context=self.task_ctx)


TaskHandler = WireTaskHandler
