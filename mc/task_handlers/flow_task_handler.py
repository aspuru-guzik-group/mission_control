from .base_task_handler import BaseTaskHandler
from mc.flow_engines import flow_engine


class FlowTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        task['data']['flow_meta'] = self.create_flow(
            task=task, task_context=task_context)

    def create_flow(self, task=None, task_context=None):
        flow = flow_engine.FlowEngine.generate_flow(
            flow_spec=task['task_params']['flow_spec'])
        serialization = flow_engine.FlowEngine.serialize_flow(flow=flow)
        flow_meta = {'serialization': serialization}
        return flow_meta

    def intermediate_tick(self, task=None, task_context=None):
        flow = self.get_flow(task=task)
        assert flow is not None
        if flow.status == 'COMPLETED':
            task['status'] = 'COMPLETED'
            task['data']['flow_data'] = flow.data
        elif flow.status == 'FAILED':
            error = str(flow.data.get('errors', '<unknown>'))
            raise Exception(error)

    def get_flow(self, task=None):
        return flow_engine.FlowEngine.deserialize_flow(
            serialized_flow=task['data']['flow_meta']['serialization'])

TaskHandler = FlowTaskHandler
