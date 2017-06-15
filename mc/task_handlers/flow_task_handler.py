from .base_task_handler import BaseTaskHandler


class FlowTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        flow = self.create_flow(task=task, task_context=task_context)
        self.tick_flow_until_has_no_pending(
            flow=flow, task=task, task_context=task_context)

    def create_flow(self, task=None, task_context=None):
        flow_engine = self.get_flow_engine(task=task, task_context=task_context)
        return flow_engine.generate_flow(
            flow_spec=task['task_params']['flow_spec'])

    def get_flow_engine(self, task=None, task_context=None):
        try: return task_context['flow_ctx']['flow_engine']
        except: raise Exception("unable to get flow_engine from task_context")

    def tick_flow_until_has_no_pending(self, flow=None, task=None,
                                       task_context=None):
        flow_engine = self.get_flow_engine(task=task, task_context=task_context)
        flow_engine.tick_flow_until_has_no_pending(
            flow=flow,
            task_context={'parent_task_context': task_context,
                          'flow_ctx': task_context.get('flow_ctx', {})}
        )
        if flow.status not in {'RUNNING', 'PENDING'}:
            self.handle_completed_flow(flow=flow, task=task,
                                       task_context=task_context)
        else:
            self.handle_incomplete_flow(flow=flow, task=task,
                                        task_context=task_context)

    def handle_completed_flow(self, flow=None, task=None, task_context=None):
        if flow.status == 'COMPLETED':
            task['status'] = 'COMPLETED'
            task['data']['flow_data'] = flow.data
        elif flow.status == 'FAILED':
            error = str(flow.data.get('errors', '<unknown>'))
            raise Exception(error)

    def handle_incomplete_flow(self, flow=None, task=None, task_context=None):
        self.persist_flow(flow=flow, task=task, task_context=task_context)

    def persist_flow(self, flow=None, task=None, task_context=None):
        flow_engine = self.get_flow_engine(task=task, task_context=task_context)
        serialization = flow_engine.serialize_flow(flow=flow)
        task['data']['flow_meta'] = {'serialization': serialization}
        task['status'] = 'RUNNING'

    def intermediate_tick(self, task=None, task_context=None):
        flow = self.get_flow(task=task, task_context=task_context)
        assert flow is not None
        self.tick_flow_until_has_no_pending(
            flow=flow, task=task, task_context=task_context)

    def get_flow(self, task=None, task_context=None):
        flow_engine = self.get_flow_engine(task=task, task_context=task_context)
        return flow_engine.deserialize_flow(
            serialized_flow=task['data']['flow_meta']['serialization'])

TaskHandler = FlowTaskHandler
