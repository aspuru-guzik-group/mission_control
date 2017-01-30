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
        create_job_fn = ctx['create_job']
        job_kwargs = {'spec': self.input['job_spec']}
        job = create_job_fn(job_kwargs=job_kwargs)
        self.data['job_uuid'] = job['uuid']
        self.status = 'RUNNING'

    def intermediate_tick(self, ctx=None):
        job = self.get_job(ctx=ctx)
        assert job is not None
        job_data = job.get('data', None) or {}
        if job['status'] == 'COMPLETED':
            self.status = 'COMPLETED'
            self.output = job_data.get('output', None)
        elif job['status'] == 'FAILED':
            self.status = 'FAILED'
            self.error = job_data.get('error', None)
        else:
            self.status = 'RUNNING'

    def get_job(self, ctx=None):
        return ctx['get_job'](uuid=self.data['job_uuid'])
