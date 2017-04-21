from mc.task_handlers.base_task_handler import BaseTaskHandler
from mc.flow_engines import flow_engine


class RunFlowTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        flow_kwargs = self.generate_flow_kwargs(
            flow_spec=task['task_params']['flow_spec'])
        created_flow = task_context['flow_ctx']['create_flow'](
            flow_kwargs=flow_kwargs)
        task['data']['flow_uuid'] = created_flow['uuid']

    def generate_flow_kwargs(self, flow_spec=None):
        flow_kwargs = {
            'label': flow_spec.get('label'),
            'serialization': flow_engine.FlowEngine.serialize_flow(
                flow=flow_engine.FlowEngine.generate_flow(flow_spec=flow_spec)
            )
        }
        return flow_kwargs

    def intermediate_tick(self, task=None, task_context=None):
        flow_ctx = task_context['flow_ctx']
        flow = self.get_flow(task=task, flow_ctx=flow_ctx)
        assert flow is not None
        if flow['status'] == 'COMPLETED':
            task['status'] = 'COMPLETED'
            task['data']['flow_data'] = flow.get('data')
        elif flow['status'] == 'FAILED':
            try: error = str(flow['serialization']['data']['errors'])
            except KeyError: error = '<unknown>'
            raise Exception(error)

    def get_flow(self, task=None, flow_ctx=None):
        return flow_ctx['get_flow'](uuid=task['data']['flow_uuid'])

TaskHandler = RunFlowTaskHandler
