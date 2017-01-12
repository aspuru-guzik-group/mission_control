import time

from flow_engines.flow_engine import FlowEngine
from flow_client.flow_client import MissionControlFlowClient as FlowClient
from job_spec_client.job_spec_client import MissionControlJobSpecClient \
        as JobClient
from .odyssey_push_job_runner import OdysseyPushJobRunner
from flow_runners.base_flow_runner import BaseFlowRunner as FlowRunner

class OdysseyPushRunner(object):
    def __init__(self, *args, run_setup=True, tick_interval=None, **kwargs):
        if run_setup: self.setup(**kwargs)
        self.tick_interval = tick_interval

    def setup(self,
              flow_generator_classes=None, 
              flow_engine=None, 
              flow_client=None, 
              job_client=None,
              job_runner=None,
              flow_runner=None,
             ):
        self.flow_generator_classes = flow_generator_classes or \
                self.generate_flow_generator_classes()
        self.flow_engine = flow_engine or self.generate_flow_engine()
        self.flow_client = flow_client or self.generate_flow_client()
        self.job_client = job_client or self.generate_job_client()
        self.job_runner = job_runner or self.generate_job_runner()
        self.flow_runner = flow_runner or self.generate_flow_runner()

    def generate_flow_generator_classes(self):
        flow_generator_classes = set()
        from a2g2.flow_generators import reaxys
        flow_generator_classes.add(reaxys.ReaxysFlowGenerator)
        return flow_generator_classes

    def generate_flow_engine(self, flow_generator_classes=None):
        flow_engine = FlowEngine()
        for flow_generator_class in (flow_generator_classes or []):
            flow_engine.register_flow_generator_class(
                flow_generator_class=flow_generator_class)
        return flow_engine

    def generate_flow_client(self): return FlowClient()

    def generate_job_client(self): return JobClient()

    def generate_job_runner(self, job_client=None): 
        return OdysseyPushJobRunner(job_client=job_client)

    def generate_flow_runner(self, flow_client=None, job_client=None,
                             flow_engine=None):
        return FlowRunner(flow_client=flow_client, job_client=job_client,
                          flow_engine=flow_engine)

    def create_flow(self, *args, flow=None, **kwargs):
        return self.flow_client.create_flow(flow=flow)

    def run(self, ntimes=None, tick_interval=None):
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep(tick_interval=tick_interval)
        else:
            while self._ticking:
                self._tick_and_sleep(tick_interval=tick_interval)

    def _tick_and_sleep(self, tick_interval=None):
        if tick_interval is None: tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        self.flow_runner.tick()
        self.job_runner.tick()
