from .base import BaseNode


class JobNode(BaseNode):
    def __init__(self, *args, create_job=None, jobs=None, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.create_job = create_job
        self.jobs = jobs

    def tick(self):
        self.increment_tick_counter()
        try: 
            if self.data['ticks'] == 1: self.initial_tick()
            else: self.intermediate_tick()
        except Exception as e:
            self.logger.exception(e)
            self.status = 'FAILED'

    def initial_tick(self):
        self.data['job_id'] = self.create_job(job_kwargs={
            'type': self.data['input']['job_type'],
            'spec': self.data['input'].get('job_spec', {})
        })
        self.status = 'RUNNING'

    def intermediate_tick(self):
        job = self.get_job()
        assert job is not None
        if job['status'] == 'COMPLETED':
            self.status = 'COMPLETED'
            self.data['outputs'] = {
                'job_data': job['data'],
            }
        elif job['status'] == 'FAILED':
            self.status = 'FAILED'
            self.data['error'] = job['data']['error']
        else:
            self.status = 'RUNNING'

    def get_job(self):
        return self.jobs[self.data['job_id']]
