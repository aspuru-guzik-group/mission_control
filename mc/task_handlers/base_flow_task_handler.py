from .base_task_handler import BaseTaskHandler


class BaseFlowTaskHandler(BaseTaskHandler):
    def get_flow_engine(self):
        try: return self.task_ctx['flow_engine']
        except Exception as exc:
            raise Exception("unable to get flow_engine") from exc

    def handle_flow_status(self, flow=None):
        if flow.status == 'COMPLETED':
            self.task['status'] = 'COMPLETED'
            self.task['data']['flow_data'] = flow.data
        elif flow.status == 'FAILED':
            error = str(flow.data.get('errors', '<unknown>'))
            raise Exception(error)
