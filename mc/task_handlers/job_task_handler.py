from .base_task_handler import BaseTaskHandler


class JobTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        task['data']['_job_meta'] = self.create_job(
            task=task, task_context=task_context)

    def intermediate_tick(self, task=None, task_context=None):
        job = self.get_job(task=task, task_context=task_context)
        assert job is not None
        if job['status'] == 'COMPLETED':
            task['data']['artifact'] = job['data'].get('artifact')
            task['data']['stdout'] = job['data'].get('stdout')
            task['status'] = 'COMPLETED'
        elif job['status'] == 'FAILED':
            try: error = job['data']['error']
            except KeyError: error = '<unknown>'
            raise Exception(error)

    def create_job(self, task=None, task_context=None):
        create_job_fn = task_context['mc.tasks.job.create_job']
        job_kwargs = {'job_spec': task['task_params']['job_spec']}
        job_meta = create_job_fn(job_kwargs=job_kwargs)
        return job_meta

    def get_job(self, task=None, task_context=None):
        get_job_fn = task_context['mc.tasks.job.get_job']
        job = get_job_fn(job_meta=task['data']['_job_meta'])
        return job

TaskHandler = JobTaskHandler
