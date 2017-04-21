from mc.task_handlers.base_task_handler import BaseTaskHandler


class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        geom = self.parse_geom()
        task['data']['geom'] = geom
        task['data']['chemthing_actions'] = self.generate_chemthing_actions(
            geom=geom, task=task)
        task['status'] = 'COMPLETED'

    def parse_geom(self):
        return {'atoms': [1,2,3]}

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
                    'props': {
                        'a2g2:prop:atoms': geom['atoms'],
                    }
                }
            }
        ]
        return chemthing_actions
