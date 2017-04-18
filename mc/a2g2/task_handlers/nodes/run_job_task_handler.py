from mc.task_handlers.base_task_handler import BaseTaskHandler


class RunJobTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        flow_ctx = task_context['flow_ctx']
        create_job_fn = flow_ctx['create_job']
        job_kwargs = {'job_spec': task['task_params']['job_spec']}
        job = create_job_fn(job_kwargs=job_kwargs)
        task['data']['job_uuid'] = job['uuid']

    def intermediate_tick(self, task=None, task_context=None):
        flow_ctx = task_context['flow_ctx']
        job = self.get_job(task=task, flow_ctx=flow_ctx)
        assert job is not None
        if job['status'] == 'COMPLETED':
            task['data']['artifact'] = job['data'].get('artifact')
            task['data']['stdout'] = job['data'].get('stdout')
            task['status'] = 'COMPLETED'
        elif job['status'] == 'FAILED':
            try: error = job['data']['error']
            except KeyError: error = '<unknown>'
            raise Exception(error)

    def get_job(self, task=None, flow_ctx=None):
        return flow_ctx['get_job'](uuid=task['data']['job_uuid'])

TaskHandler = RunJobTaskHandler
