from .base_task_handler import BaseTaskHandler


class JobTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        task['data']['_job_meta'] = self.create_job(
            task=task, task_context=task_context)

    def intermediate_tick(self, task=None, task_context=None):
        job = self.get_job(task=task, task_context=task_context)
        assert job is not None
        job_data = job.get('data', {})
        if job['status'] == 'COMPLETED':
            task['data']['artifact'] = job_data.get('artifact')
            task['data']['std_logs'] = job_data.get('std_logs')
            task['status'] = 'COMPLETED'
        elif job['status'] == 'FAILED':
            error = job_data.get('error', '<unknown job error>')
            raise Exception(error)

    def create_job(self, task=None, task_context=None):
        create_job_fn = task_context['mc.tasks.job.create_job']
        job_kwargs = {'job_spec': task['task_params'].get('job_spec')}
        job_meta = create_job_fn(job_kwargs=job_kwargs)
        return job_meta

    def get_job(self, task=None, task_context=None):
        get_job_fn = task_context['mc.tasks.job.get_job']
        job = get_job_fn(job_meta=task['data']['_job_meta'])
        return job

TaskHandler = JobTaskHandler
