from mc.task_handlers.base_task_handler import BaseTaskHandler

from .. import utils as _a2g2_db_utils

class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        a2g2_client = _a2g2_db_utils.generate_a2g2_client(
            cfg=task_context['flow_ctx']['cfg'])
        a2g2_client.post_chemthing_actions(
            chemthing_actions=task['task_params']['chemthing_actions'])
        task['status'] = 'COMPLETED'
