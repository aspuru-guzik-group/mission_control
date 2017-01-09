from .base import BaseTask


class FlowTask(BaseTask):
    def tick(self, *args, ctx=None, **kwargs):
        self.increment_tick_counter()
        try: 
            if self.data['ticks'] == 1: self.initial_tick(ctx=ctx)
            else: self.intermediate_tick(ctx=ctx)
        except Exception as e:
            self.mark_as_failed(error=e)
            raise e

    def initial_tick(self, ctx=None):
        create_flow = ctx['create_flow']
        self.data['flow_uuid'] = create_flow(flow_kwargs={
            'flow_spec': self.input.get('flow_spec', {})
        })
        self.status = 'RUNNING'

    def intermediate_tick(self, ctx=None):
        flow = self.get_flow(ctx=ctx)
        assert flow is not None
        if flow['status'] == 'COMPLETED':
            self.status = 'COMPLETED'
            self.output = flow.get('output', None)
        elif flow['status'] == 'FAILED':
            self.status = 'FAILED'
            self.error = flow.get('error', None)
        else:
            self.status = 'RUNNING'

    def get_flow(self, ctx=None):
        return ctx['get_flow'](flow_uuid=self.data['flow_uuid'])
