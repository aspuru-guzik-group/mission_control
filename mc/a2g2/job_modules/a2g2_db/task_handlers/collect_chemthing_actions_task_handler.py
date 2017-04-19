from mc.task_handlers.base_task_handler import BaseTaskHandler


class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        task['data']['chemthing_actions'] = self.collect_chemthing_actions(
            task_context=task_context)
        task['status'] = 'COMPLETED'

    def collect_chemthing_actions(self, task_context=None):
        chemthing_actions = []
        for node_task in task_context['node']['node_tasks']:
            chemthing_actions.extend(node_task.get('data', {}).get(
                'chemthing_actions', []))
        return chemthing_actions
