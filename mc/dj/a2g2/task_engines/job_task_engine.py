from .base_task_engine import BaseTaskEngine


class JobTaskEngine(BaseTaskEngine):
    def tick_task(self, *args, task=None, ctx=None, **kwargs):
        self.increment_task_tick_counter(task=task)
        try: 
            tick_kwargs = {'task': task, 'ctx': ctx}
            if task['data']['ticks'] == 1: self.initial_tick(**tick_kwargs)
            else: self.intermediate_tick(**tick_kwargs)
        except Exception as e:
            self.mark_task_as_failed(task=task, error=e)
            raise e

    def initial_tick(self, task=None, ctx=None):
        create_job_fn = ctx['create_job']
        job_kwargs = {'spec': task['input']['job_spec']}
        job = create_job_fn(job_kwargs=job_kwargs)
        task['data']['job_uuid'] = job['uuid']
        task['status'] = 'RUNNING'

    def intermediate_tick(self, task=None, ctx=None):
        job = self.get_job(task=task, ctx=ctx)
        assert job is not None
        job_data = job.get('data', None) or {}
        if job['status'] == 'COMPLETED':
            task['status'] = 'COMPLETED'
            task['output'] = job_data.get('output', None)
        elif job['status'] == 'FAILED':
            self.mark_task_as_failed(task=task,
                                     error=job_data.get('error', None))
        else:
            task['status'] = 'RUNNING'

    def get_job(self, task=None, ctx=None):
        return ctx['get_job'](uuid=task['data']['job_uuid'])
