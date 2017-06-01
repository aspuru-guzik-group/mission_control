from mc.task_handlers.base_task_handler import BaseTaskHandler


class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        task['data']['chemthing_actions'] = self.generate_chemthing_actions(
            task=task)
        task['status'] = 'COMPLETED'

    def generate_chemthing_actions(self, geom=None, task=None):
        chemthing_actions = [
            {
                'key': task['task_params']['key'],
                'updates': {
                    'tags': {
                        task['task_params']['key']: True
                    },
                    'types': {
                        'a2g2:type:mol3d': True,
                    },
                    'props': task['task_params']['props']
                }
            }
        ]
        return chemthing_actions
