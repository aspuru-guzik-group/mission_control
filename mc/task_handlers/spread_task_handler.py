from mc.utils import context_value_setter

from .base_proxying_task_handler import BaseProxyingTaskHandler


class SpreadTaskHandler(BaseProxyingTaskHandler):
    DEFAULT_CONTAINER_TASK_TYPE = 'mc.tasks.flow'

    def validate_task_params(self):
        assert self.get_items()
        assert self.get_mapping_params()

    def get_items(self):
        return self.task['task_params']['items']

    def get_mapping_params(self):
        return self.task['task_params']['mapping_params']

    def get_container_task_type(self):
        return self.task['task_params'].get(
            'container_task_type', self.DEFAULT_CONTAINER_TASK_TYPE)

    def generate_proxied_task(self):
        return {
            'task_type': self.get_container_task_type(),
            'task_params': {
                'flow_spec': self.generate_flow_spec()
            }
        }

    def generate_flow_spec(self):
        flow_spec = {
            'label': self.get_flow_label(),
            'tasks': self.map_items_to_tasks()
        }
        return flow_spec

    def get_flow_label(self): return self.task['key'] + ':spread_subflow'

    def map_items_to_tasks(self):
        return context_value_setter.execute_mapping_transform(
            items=self.get_items(),
            transform={'type': 'mapping', 'params': self.get_mapping_params()},
            context=self.task_ctx
        )

    def on_proxied_task_finished(self):
        super().on_proxied_task_finished()
        self.task['data'] = self.proxied_task.get('flow_data')


TaskHandler = SpreadTaskHandler
