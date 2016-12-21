import logging


class ConfgenNode(object):
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
        assert self.data['inputs']['smiles'] is not None
        self.data['job_id'] = self.create_job(job_kwargs={
            'type': 'confgen',
            'spec': {
                'smiles': self.data['inputs']['smiles']
            }
        })
        self.status = 'RUNNING'

    def intermediate_tick(self):
        job = self.get_job()
        assert job is not None
        if job['status'] == 'COMPLETED':
            self.status = 'COMPLETED'
            self.data['outputs'] = {
                'conformer_dir_uri': job['data']['output_dir_uri']
            }
        elif job['status'] == 'FAILED':
            self.status = 'FAILED'
            self.data['error'] = job['data']['error']
        else:
            self.status = 'RUNNING'

    def get_job(self):
        return self.jobs[self.data['job_id']]
