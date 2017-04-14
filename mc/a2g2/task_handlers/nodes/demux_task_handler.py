from copy import deepcopy

from mc.mc_utils.context_value_setter import set_context_values
from mc.task_handlers.base_task_handler import BaseTaskHandler


class DemuxTaskHandler(BaseTaskHandler):
    ROOT_NODE_KEY = 'ROOT'

    def initial_tick(self, task=None, task_context=None):
        flow_task = self.generate_flow_task(task=task,
                                            task_context=task_context)
        task_context['node']['node_tasks'].append(flow_task)
        task['status'] = 'COMPLETED'

    def generate_flow_task(self, task=None, task_context=None):
        label = task['task_key'] + '__flow'
        flow_task = {
            'task_key': label,
            'task_type': 'a2g2.tasks.nodes.run_flow',
            'task_params': {
                'flow_spec': self.generate_flow_spec(task=task,
                                                     task_context=task_context,
                                                     label=label)
            }
        }
        return flow_task

    def generate_flow_spec(self, task=None, task_context=None, label=None):
        flow_spec = {}
        flow_spec['label'] = label
        demux_flow_params = task['task_params'].get('demux_flow_params', {})
        flow_spec['data'] = demux_flow_params.get('data', {})
        flow_spec['cfg'] = demux_flow_params.get('cfg', {})
        flow_spec['node_specs'] = [
            self.generate_node_spec(item=item, index=index,
                                    task=task, task_context=task_context)
            for index, item in enumerate(task['task_params'].get('items', []))
        ]
        return flow_spec

    def generate_node_spec(self, item=None, index=None, task=None,
                           task_context=None):
        node_spec = deepcopy(task['task_params']['node_spec_template'])
        node_spec.setdefault('node', {})
        node_spec['node']['node_key'] = "{task_key}_{index}".format(
            task_key=task['task_key'], index=index)
        node_spec_context = {'demux_item': item, 'demux_index': index,
                             'node_spec': node_spec}
        set_context_values(value_specs=task['task_params'].get('substitutions'),
                           context=node_spec_context)
        node_spec['precursor_keys'] = [self.ROOT_NODE_KEY]
        return node_spec
