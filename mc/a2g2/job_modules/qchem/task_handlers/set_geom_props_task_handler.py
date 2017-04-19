from mc.task_handlers.base_task_handler import BaseTaskHandler


class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        task['data']['chemthing_actions'] = self.generate_chemthing_actions(
            task=task)
        task['status'] = 'COMPLETED'

    def generate_chemthing_actions(self, geom=None, task=None):
        parent_computation = task['task_params']['parent_computation']
        key = self.generate_chemthing_key(parent_computation=parent_computation)
        chemthing_actions = [
            {
                'key': key,
                'updates': {
                    'keys': {key: True},
                    'types': {
                        'a2g2:type:mol3d': True,
                    },
                    'props': task['task_params']['props']
                }
            }
        ]
        return chemthing_actions

    def generate_chemthing_key(self, parent_computation=None):
        return "{parent_namespace}:optimized_geometry".format(
            parent_namespace=parent_computation['uuid'])
