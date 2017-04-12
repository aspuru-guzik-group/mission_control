from copy import deepcopy

from mc.mc_utils.context_value_setter import set_context_values
from mc.task_handlers.base_task_handler import BaseTaskHandler


class DemuxTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        flow_task = self.generate_flow_task(task=task,
                                            task_context=task_context)
        task_context['node']['node_tasks'].append(flow_task)
        task['status'] = 'COMPLETED'

    def generate_flow_task(self, task=None, task_context=None):
        flow_task = {
            'task_type': 'a2g2:task:run_flow',
            'task_params': {
                'flow_spec': self.generate_flow_spec(task=task,
                                                     task_context=task_context)
            }
        }
        return flow_task

    def generate_flow_spec(self, task=None, task_context=None):
        flow_spec = {}
        flow_spec['node_specs'] = [
            self.generate_node_spec(item=item, index=index,
                                    task=task, task_context=task_context)
            for index, item in enumerate(task['task_params'].get('items', []))
        ]
        return flow_spec

    def generate_node_spec(self, item=None, index=None, task=None,
                           task_context=None):
        node_spec = deepcopy(task['task_params']['node_spec_template_obj'])
        node_spec_context = {'item': item, 'index': index,
                             'node_spec': node_spec}
        set_context_values(value_specs=task['task_params'].get('substitutions'),
                           context=node_spec_context)
        return node_spec
