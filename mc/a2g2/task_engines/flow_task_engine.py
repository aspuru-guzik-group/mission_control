import json
from .base_task_engine import BaseTaskEngine


class FlowTaskEngine(BaseTaskEngine):
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
        flow = ctx['create_flow'](flow={
            'spec': json.dumps(task['input'].get('flow_spec', {}))
        })
        task['data']['flow_uuid'] = flow['uuid']
        task['status'] = 'RUNNING'

    def intermediate_tick(self, task=None, ctx=None):
        flow = self.get_flow(task=task, ctx=ctx)
        assert flow is not None
        if flow['status'] == 'COMPLETED':
            task['status'] = 'COMPLETED'
            task['output'] = flow.get('output', None)
        elif flow['status'] == 'FAILED':
            self.mark_task_as_failed(task=task, error=flow.get('error', None))
        else:
            task['status'] = 'RUNNING'

    def get_flow(self, task=None, ctx=None):
        return ctx['get_flow'](uuid=task['data']['flow_uuid'])
