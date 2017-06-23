from mc.task_handlers.base_flow_task_handler import BaseFlowTaskHandler


class FlowTaskHandler(BaseFlowTaskHandler):
    def initial_tick(self):
        flow_meta = self.create_flow_record()
        self.task['data']['_flow_task_flow_meta'] = flow_meta

    def create_flow_record(self):
        flow = self.create_flow()
        create_record_fn = self.task_ctx['mc.tasks.flow.create_flow_record']
        flow_engine = self.get_flow_engine()
        flow_record_kwargs = flow_engine.flow_to_flow_dict(flow=flow)
        flow_record = create_record_fn(flow_kwargs=flow_record_kwargs)
        flow_meta = {'key': flow_record['key']}
        return flow_meta

    def create_flow(self):
        flow_engine = self.get_flow_engine()
        flow = flow_engine.flow_spec_to_flow(
            flow_spec=self.task['task_params']['flow_spec'])
        flow.depth = self.task_ctx['flow'].depth + 1
        flow.data['parent_key'] = self.task_ctx['flow'].key
        return flow

    def intermediate_tick(self):
        flow = self.get_flow()
        self.handle_flow_status(flow=flow)

    def get_flow(self):
        get_flow_fn = self.task_ctx['mc.tasks.flow.get_flow_record']
        flow_record = get_flow_fn(
            flow_meta=self.task['data']['_flow_task_flow_meta'])
        flow_engine = self.get_flow_engine()
        return flow_engine.flow_dict_to_flow(flow_dict=flow_record)

TaskHandler = FlowTaskHandler
