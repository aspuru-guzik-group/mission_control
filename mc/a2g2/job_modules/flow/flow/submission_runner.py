from mc.flow_engines import flow_engine

from mc.a2g2.task_handlers.a2g2_task_handler import A2G2TaskHandler
from ...a2g2_common.base_submission_runner import BaseSubmissionRunner


class SubmissionRunner(BaseSubmissionRunner):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.flow_engine = self.generate_flow_engine()

    def generate_flow_engine(self):
        return flow_engine.FlowEngine(task_handler=A2G2TaskHandler())

    def run_submission(self):
        self.run_flow(flow=self.generate_flow(),
                      flow_ctx=self.generate_flow_ctx())

    def generate_flow(self):
        return self.flow_engine.generate_flow(flow_spec=self.get_flow_spec())

    def get_flow_spec(self):
        return self.job['job_spec']['job_params']['flow_spec']

    def generate_flow_ctx(self):
        return {
            'job': self.job,
            'cfg': self.cfg,
            'submission': self.submission
        }

    def run_flow(self, flow=None, flow_ctx=None):
        self.flow_engine.run_flow(flow=flow, flow_ctx=flow_ctx)
        if flow.status == 'FAILED':
            raise Exception("Flow failed. flow.data['errors']: {}".format(
                flow.data.get('errors')))

def run_job_submission(*args, job=None, cfg=None, submission=None, **kwargs):
    runner = SubmissionRunner(job=job, cfg=cfg, submission=submission)
    return runner.run_submission()
