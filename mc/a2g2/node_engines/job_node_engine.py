from .base_node_engine import BaseNodeEngine


class JobNodeEngine(BaseNodeEngine):
    def tick_node(self, *args, node=None, ctx=None, **kwargs):
        self.increment_node_tick_counter(node=node)
        try: 
            tick_kwargs = {'node': node, 'ctx': ctx}
            if node['data']['ticks'] == 1: self.initial_tick(**tick_kwargs)
            else: self.intermediate_tick(**tick_kwargs)
        except Exception as e:
            self.mark_node_as_failed(node=node, error=e)
            raise e

    def initial_tick(self, node=None, ctx=None):
        create_job_fn = ctx['create_job']
        job_kwargs = {'job_spec': node['input']['job_spec']}
        job = create_job_fn(job_kwargs=job_kwargs)
        node['data']['job_uuid'] = job['uuid']
        node['status'] = 'RUNNING'

    def intermediate_tick(self, node=None, ctx=None):
        job = self.get_job(node=node, ctx=ctx)
        assert job is not None
        job_data = job.get('data', None) or {}
        if job['status'] == 'COMPLETED':
            node['status'] = 'COMPLETED'
            node['output'] = job_data.get('output', None)
        elif job['status'] == 'FAILED':
            self.mark_node_as_failed(node=node,
                                     error=job_data.get('error', None))
        else:
            node['status'] = 'RUNNING'

    def get_job(self, node=None, ctx=None):
        return ctx['get_job'](uuid=node['data']['job_uuid'])
