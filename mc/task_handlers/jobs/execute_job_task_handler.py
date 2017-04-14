from ..base_task_handler import BaseTaskHandler


class ExecuteJobTaskHandler(BaseTaskHandler):
    def __init__(self, *args, execution_client=None, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.execution_client = execution_client

    def initial_tick(self, task=None, task_context=None):
        execution_meta = self.execution_client.start_execution(
            submission=task['task_params']['submission'])
        task['data']['execution_meta'] = execution_meta

    def intermediate_tick(self, task=None, task_context=None):
        execution_state = self.execution_client.get_execution_state(
            execution_meta=task['data']['execution_meta'])
        task['data']['execution_state'] = execution_state
        self.handle_execution_state(execution_state=execution_state, task=task,
                                    job=task_context['job'])

    def handle_execution_state(self, execution_state=None, task=None, job=None):
        subhandler_kwargs = {'execution_state': execution_state, 'task': task,
                             'job': job}
        run_status = execution_state['run_status']
        if run_status == 'COMPLETED':
            self.handle_completed_execution_state(**subhandler_kwargs)
        elif run_status == 'FAILED':
            self.handle_failed_execution_state(**subhandler_kwargs)

    def handle_completed_execution_state(self, execution_state=None, task=None,
                                         job=None):
        for attr in ['artifact', 'stdout']:
            task['data'][attr] = execution_state.get(attr)
        task['status'] = 'COMPLETED'

    def handle_failed_execution_state(self, execution_state=None, task=None,
                                      job=None):
        task['status'] = 'FAILED'
        task['data']['error'] = execution_state.get(
            'error', '<unknown execution error>')
