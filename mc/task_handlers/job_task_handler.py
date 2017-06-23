from .base_task_handler import BaseTaskHandler


class JobTaskHandler(BaseTaskHandler):
    def initial_tick(self):
        print("itick")
        self.task['data']['_job_task_job_meta'] = self.create_job()

    def intermediate_tick(self):
        job = self.get_job()
        job_data = job.get('data', {})
        if job['status'] == 'COMPLETED':
            self.task['data']['artifact'] = job_data.get('artifact')
            self.task['data']['std_logs'] = job_data.get('std_logs')
            self.task['status'] = 'COMPLETED'
        elif job['status'] == 'FAILED':
            error = job_data.get('error', '<unknown job error>')
            raise Exception(error)

    def create_job(self):
        create_job_fn = self.task_ctx['mc.tasks.job.create_job']
        job_kwargs = {
            'job_spec': self.task['task_params'].get('job_spec'),
            'data': {'parent_key': self.task_ctx['flow'].key}
        }
        job_meta = create_job_fn(job_kwargs=job_kwargs)
        return job_meta

    def get_job(self):
        get_job_fn = self.task_ctx['mc.tasks.job.get_job']
        job = get_job_fn(job_meta=self.task['data']['_job_task_job_meta'])
        return job

TaskHandler = JobTaskHandler
