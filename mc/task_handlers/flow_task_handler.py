from mc.task_handlers.base_flow_task_handler import BaseFlowTaskHandler


class FlowTaskHandler(BaseFlowTaskHandler):
    def initial_tick(self):
        flow_meta = self.create_flow_record()
        self.task['data']['_flow_task_flow_meta'] = flow_meta

    def create_flow_record(self):
        flow = self.create_flow()
        flow_record = self.flow_record_client.create_flow_record(
            flow_kwargs=self.flow_engine.flow_to_flow_dict(flow=flow))
        flow_meta = {'key': flow_record['key']}
        return flow_meta

    def create_flow(self):
        flow = self.flow_engine.flow_spec_to_flow(
            flow_spec=self.task['task_params']['flow_spec'])
        flow.depth = self.task_ctx['flow'].depth + 1
        flow.data['parent_key'] = self.task_ctx['flow'].key
        return flow

    @property
    def flow_record_client(self):
        if not hasattr(self, '_flow_record_client'):
            self._flow_record_client = self.task_ctx['mc.flow_record_client']
        return self._flow_record_client

    @flow_record_client.setter
    def flow_record_client(self, value): self._flow_record_client = value

    def intermediate_tick(self):
        flow = self.get_flow()
        self.handle_flow_status(flow=flow)

    def get_flow(self):
        flow_record = self.get_flow_record()
        return self.flow_engine.flow_dict_to_flow(flow_dict=flow_record)

    def get_flow_record(self):
        return self.flow_record_client.get_flow_record(
            flow_meta=self.task['data']['_flow_task_flow_meta'])

TaskHandler = FlowTaskHandler
