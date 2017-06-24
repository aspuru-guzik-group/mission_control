from .base_task_handler import BaseTaskHandler


class BaseFlowTaskHandler(BaseTaskHandler):
    class FlowTaskError(Exception): pass

    @property
    def flow_engine(self):
        if not hasattr(self, '_flow_engine'):
            try: self._flow_engine = self.task_ctx['flow_engine']
            except Exception as exc:
                raise self.FlowTaskError("unable to get flow_engine") from exc
        return self._flow_engine

    @flow_engine.setter
    def flow_engine(self, value): self._flow_engine = value

    def handle_flow_status(self, flow=None):
        try:
            if flow.status == 'COMPLETED':
                self.task['status'] = 'COMPLETED'
                self.task['data']['flow_data'] = flow.data
            elif flow.status == 'FAILED':
                error = str(flow.data.get('errors', '<unknown>'))
                raise self.FlowTaskError(error)
        except: raise
        finally: 
            if flow.status in {'COMPLETED', 'FAILED'}:
                self.on_flow_finished(flow=flow)

    def on_flow_finished(self, flow=None): pass
