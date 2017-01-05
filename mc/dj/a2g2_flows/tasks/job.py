from .base import BaseTask


class JobTask(BaseTask):
    def __init__(self, *args, **kwargs): super().__init__(self, *args, **kwargs)

    def tick(self, *args, ctx=None, **kwargs):
        self.increment_tick_counter()
        try: 
            if self.data['ticks'] == 1: self.initial_tick(ctx=ctx)
            else: self.intermediate_tick(ctx=ctx)
        except Exception as e:
            self.logger.exception(e)
            self.status = 'FAILED'

    def initial_tick(self, ctx=None):
        create_job = ctx['create_job']
        self.data['job_id'] = create_job(job_kwargs={
            'type': self.input['job_type'],
            'spec': self.input.get('job_spec', {})
        })
        self.status = 'RUNNING'

    def intermediate_tick(self, ctx=None):
        job = self.get_job(ctx=ctx)
        assert job is not None
        if job['status'] == 'COMPLETED':
            self.status = 'COMPLETED'
            self.output = job.get('output', None)
        elif job['status'] == 'FAILED':
            self.status = 'FAILED'
            self.error = job.get('error', None)
        else:
            self.status = 'RUNNING'

    def get_job(self, ctx=None):
        return ctx['get_job'](job_id=self.data['job_id'])
