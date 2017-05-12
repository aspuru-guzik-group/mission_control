from mc.flow_engines.flow_engine import FlowEngine
from .base_flow_runner import BaseFlowRunner

from a2g2.task_handlers import a2g2_task_handler


class FlowRunner(object):
    def __init__(self, mc_client=None, flow_client=None, task_handler=None,
                 flow_ctx=None):
        self.mc_client = mc_client
        self.flow_client = flow_client or self.mc_client
        self.task_handler = task_handler or a2g2_task_handler.A2G2TaskHandler
        self.flow_ctx = self._decorate_flow_ctx(flow_ctx=flow_ctx)

        self._base_flow_runner = self._generate_base_flow_runner()
        self._tick_counter = 0

    def _generate_base_flow_runner(self, flow_runner_kwargs=None):
        return BaseFlowRunner(flow_client=self.flow_client,
                              flow_engine=self._generate_flow_engine(),
                              flow_ctx=self.flow_ctx)

    def _generate_flow_engine(self):
        return FlowEngine(task_handler=self.task_handler, max_msg_len=4048)

    def _decorate_flow_ctx(self, flow_ctx=None):
        flow_ctx = flow_ctx or {}
        decorated_flow_ctx = {
            'create_job': self.mc_client.create_job,
            'get_job': self.mc_client.get_job_by_uuid,
            'create_flow': self.mc_client.create_flow,
            'get_flow': self.mc_client.get_flow_by_uuid,
            **flow_ctx
        }
        return decorated_flow_ctx

    def tick(self):
        self._tick_counter += 1
        self._base_flow_runner.tick()
