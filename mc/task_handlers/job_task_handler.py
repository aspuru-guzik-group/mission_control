from .base_task_handler import BaseTaskHandler


class JobTaskHandler(BaseTaskHandler):
    def initial_tick(self):
        self.task['data']['_job_task_job_meta'] = self.create_job_record()

    def intermediate_tick(self):
        job_record = self.get_job_record()
        job_data = job_record.get('data', {})
        job_status = job_record.get('status')
        if job_status == 'COMPLETED':
            self.task['data']['artifact'] = job_data.get('artifact')
            self.task['data']['std_logs'] = job_data.get('std_logs')
            self.task['status'] = 'COMPLETED'
        elif job_status == 'FAILED':
            error = job_data.get('error', '<unknown job error>')
            raise Exception(error)

    def create_job_record(self):
        job_kwargs = {
            'job_spec': self.task['task_params'].get('job_spec'),
            'data': {'parent_key': self.task_ctx['flow'].key}
        }
        job_meta = self.job_record_client.create_job_record(
            job_kwargs=job_kwargs)
        return job_meta

    @property
    def job_record_client(self):
        if not hasattr(self, '_job_record_client'):
            self._job_record_client = self.task_ctx['mc.job_record_client']
        return self._job_record_client

    @job_record_client.setter
    def job_record_client(self, value): self._job_record_client = value

    def get_job_record(self):
        job_record = self.job_record_client.get_job_record(
            job_meta=self.task['data']['_job_task_job_meta'])
        return job_record

TaskHandler = JobTaskHandler
