import logging


class JobNode(object):
    def __init__(self, *args, create_job=None, status=None, data=None,
                 jobs=None, logger=None, **kwargs):
        self.create_job = create_job
        self.status = status
        self.data = data
        self.jobs = jobs
        self.logger = logger or logging

    def tick(self):
        self.increment_tick_counter()
        try: 
            if self.data['ticks'] == 1: self.initial_tick()
            else: self.intermediate_tick()
        except Exception as e:
            self.logger.exception(e)
            self.status = 'FAILED'

    def increment_tick_counter(self):
        if 'ticks' not in self.data: self.data['ticks'] = 0
        self.data['ticks'] += 1

    def initial_tick(self):
        self.data['job_id'] = self.create_job(job_kwargs={
            'type': self.data['inputs']['job_type'],
            'spec': self.data['inputs'].get('job_spec', {})
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
