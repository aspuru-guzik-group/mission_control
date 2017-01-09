from .base import BaseTask


class JobTask(BaseTask):
    def tick(self, *args, ctx=None, **kwargs):
        self.increment_tick_counter()
        try: 
            if self.data['ticks'] == 1: self.initial_tick(ctx=ctx)
            else: self.intermediate_tick(ctx=ctx)
        except Exception as e:
            self.mark_as_failed(error=e)
            raise e

    def initial_tick(self, ctx=None):
        create_job = ctx['create_job']
        job_kwargs = {'job_spec': self.input['job_spec']}
        self.data['job_id'] = create_job(job_kwargs=job_kwargs)
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
