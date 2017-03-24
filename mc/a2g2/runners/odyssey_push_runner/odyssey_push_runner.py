import time

import requests

from mc.mc_client.mission_control_client import MissionControlClient
from mc.flow_engines.flow_engine import FlowEngine
from mc.flow_runners.base_flow_runner import BaseFlowRunner as FlowRunner
from mc.job_runners.base_job_runner import BaseJobRunner as JobRunner
from .odyssey_job_submission_factory import OdysseyJobSubmissionFactory as \
        JobSubmissionFactory


class OdysseyPushRunner(object):
    def __init__(self, *args, request_client=None, run_setup=True,
                 tick_interval=None, mc_server_url=None, **setup_kwargs):
        self.request_client = request_client or requests
        self.tick_interval = tick_interval
        self.mc_server_url = mc_server_url
        if run_setup: self.setup(**setup_kwargs)
        self.tick_counter = 0

    def setup(self,
              task_handler=None,
              flow_generator_classes=None, 
              flow_engine=None, 
              mc_client=None, 
              job_submission_factory=None,
              job_runner=None,
              job_runner_kwargs=None,
              tick_ctx=None,
              flow_runner=None,
              **kwargs
             ):
        self.task_handler = task_handler or self.generate_task_handler()
        self.flow_generator_classes = flow_generator_classes or \
                self.generate_flow_generator_classes()
        self.flow_engine = flow_engine or self.generate_flow_engine()
        self.mc_client = mc_client or self.generate_mc_client()
        self.job_submission_factory = job_submission_factory or \
                self.generate_job_submission_factory()
        self.job_runner = job_runner or self.generate_job_runner(
            job_runner_kwargs=job_runner_kwargs)
        self.tick_ctx = self.decorate_tick_ctx(tick_ctx=tick_ctx)
        self.flow_runner = flow_runner or self.generate_flow_runner()

    def generate_task_handler(self): pass

    def generate_flow_generator_classes(self):
        flow_generator_classes = set()
        return flow_generator_classes

    def generate_flow_engine(self):
        flow_engine = FlowEngine(task_handler=self.task_handler)
        for flow_generator_class in (self.flow_generator_classes or []):
            flow_engine.register_flow_generator_class(
                flow_generator_class=flow_generator_class)
        return flow_engine

    def generate_mc_client(self):
        return MissionControlClient(base_url=self.mc_server_url,
                                    request_client=self.request_client)

    def generate_job_submission_factory(self):
        return JobSubmissionFactory()

    def generate_job_runner(self, job_runner_kwargs=None): 
        return JobRunner(
            task_handler=self.task_handler,
            job_submission_factory=self.job_submission_factory,
            job_client=self.mc_client,
            **(job_runner_kwargs or {})
        )

    def decorate_tick_ctx(self, tick_ctx=None):
        tick_ctx = tick_ctx or {}
        decorated_tick_ctx = {
            'create_job': self.mc_client.create_job,
            'get_job': self.mc_client.fetch_job_by_uuid,
            'create_flow': self.mc_client.create_flow,
            'get_flow': self.mc_client.fetch_flow_by_uuid,
            **tick_ctx
        }
        return decorated_tick_ctx

    def generate_flow_runner(self):
        return FlowRunner(flow_client=self.mc_client,
                          flow_engine=self.flow_engine,
                          tick_ctx=self.tick_ctx)

    def create_flow_record(self, *args, flow_record=None, **kwargs):
        return self.mc_client.create_flow(flow=flow_record)

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
        self.tick_counter += 1
        self.flow_runner.tick()
        self.job_runner.tick()
