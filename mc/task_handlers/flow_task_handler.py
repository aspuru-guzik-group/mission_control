from mc.task_handlers.base_task_handler import BaseTaskHandler


class FlowTaskHandler(BaseTaskHandler):
    def initial_tick(self):
        flow_meta = self.create_flow()
        self.task['data']['flow_meta'] = flow_meta

    def create_flow(self):
        flow_spec = self.task['task_params']['flow_spec']
        flow_record_kwargs = self.generate_flow_record_kwargs(
            flow_spec=flow_spec)
        flow_record = self.task_ctx['flow_ctx']['create_flow_record'](
            kwargs=flow_record_kwargs)
        flow_meta = {'uuid': flow_record['uuid']}
        return flow_meta

    def generate_flow_record_kwargs(self, flow_spec=None):
        flow_engine = self.get_flow_engine()
        flow_record_kwargs = {
            'label': flow_spec.get('label'),
            'serialization': flow_engine.FlowEngine.serialize_flow(
                flow=flow_engine.FlowEngine.generate_flow(flow_spec=flow_spec)
            )
        }
        return flow_record_kwargs

    def intermediate_tick(self):
        flow = self.get_flow()
        if flow.status == 'COMPLETED':
            self.task['status'] = 'COMPLETED'
            self.task['data']['flow_data'] = flow.data
        elif flow.status == 'FAILED':
            error = str(flow.data.get('errors', '<unknown>'))
            raise Exception(error)

    def get_flow(self):
        flow_ctx = self.task_context['flow_ctx']
        flow_record = flow_ctx['get_flow'](
            uuid=self.task['data']['flow_meta']['uuid'])
        flow_engine = self.get_flow_engine()
        return flow_engine.deserialize_flow(
            serialized_flow=flow_record.get('serialization'))

TaskHandler = FlowTaskHandler
