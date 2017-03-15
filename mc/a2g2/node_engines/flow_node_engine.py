import json
from .base_node_engine import BaseNodeEngine


class FlowNodeEngine(BaseNodeEngine):
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
        flow = ctx['create_flow'](flow={
            'spec': json.dumps(node['input'].get('flow_spec', {}))
        })
        node['data']['flow_uuid'] = flow['uuid']
        node['status'] = 'RUNNING'

    def intermediate_tick(self, node=None, ctx=None):
        flow = self.get_flow(node=node, ctx=ctx)
        assert flow is not None
        if flow['status'] == 'COMPLETED':
            node['status'] = 'COMPLETED'
            node['output'] = flow.get('output', None)
        elif flow['status'] == 'FAILED':
            self.mark_node_as_failed(node=node, error=flow.get('error', None))
        else:
            node['status'] = 'RUNNING'

    def get_flow(self, node=None, ctx=None):
        return ctx['get_flow'](uuid=node['data']['flow_uuid'])
