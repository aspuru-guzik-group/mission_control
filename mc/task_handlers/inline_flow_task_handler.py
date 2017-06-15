from .base_task_handler import BaseTaskHandler


class InlineFlowTaskHandler(BaseTaskHandler):
    def initial_tick(self):
        self.tick_flow_until_has_no_pending(flow=self.create_flow())

    def create_flow(self):
        return self.get_flow_engine().generate_flow(
            flow_spec=self.task['task_params']['flow_spec'])

    def get_flow_engine(self):
        try: return self.task_ctx['flow_engine']
        except Exception as exc:
            raise Exception("unable to get flow_engine") from exc

    def tick_flow_until_has_no_pending(self, flow=None):
        flow_engine = self.get_flow_engine()
        flow_engine.tick_flow_until_has_no_pending(
            flow=flow,
            task_ctx={'parent_task_ctx': self.task_ctx,
                      **{k: v for k, v in self.task_ctx.items()
                         if k not in {'flow', 'task'}}
                     }
        )
        self.persist_flow(flow=flow)
        if flow.status not in {'RUNNING', 'PENDING'}:
            self.handle_completed_flow(flow=flow)

    def handle_completed_flow(self, flow=None):
        if flow.status == 'COMPLETED':
            self.task['status'] = 'COMPLETED'
            self.task['data']['flow_data'] = flow.data
        elif flow.status == 'FAILED':
            error = str(flow.data.get('errors', '<unknown>'))
            raise Exception(error)

    def persist_flow(self, flow=None):
        flow_engine = self.get_flow_engine()
        serialization = flow_engine.serialize_flow(flow=flow)
        self.task['data']['_flow_task_flow_meta'] = \
            {'serialization': serialization}
        self.task['status'] = 'RUNNING'

    def intermediate_tick(self):
        self.tick_flow_until_has_no_pending(flow=self.get_flow())

    def get_flow(self):
        flow_engine = self.get_flow_engine()
        return flow_engine.deserialize_flow(
            serialized_flow=(self.task['data']['_flow_task_flow_meta']
                             ['serialization'])
        )

TaskHandler = InlineFlowTaskHandler
