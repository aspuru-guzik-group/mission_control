import time

import requests

from flow_engines.flow_engine import FlowEngine
from flow_client.flow_client import MissionControlFlowClient as FlowClient
from job_client.job_client import MissionControlJobClient \
        as JobClient
from .odyssey_push_job_runner import OdysseyPushJobRunner
from flow_runners.base_flow_runner import BaseFlowRunner as FlowRunner
from .odyssey_job_dir_factory import OdysseyJobDirFactory as JobDirFactory


class OdysseyPushRunner(object):
    def __init__(self, *args, request_client=None, run_setup=True,
                 tick_interval=None, odyssey_user=None, odyssey_host=None,
                 job_server_url=None, flow_server_url=None, **setup_kwargs):
        self.request_client = request_client or requests
        self.tick_interval = tick_interval
        self.odyssey_user = odyssey_user
        self.odyssey_host = odyssey_host
        self.job_server_url = job_server_url
        self.flow_server_url = flow_server_url
        if run_setup: self.setup(**setup_kwargs)

    def setup(self,
              flow_generator_classes=None, 
              flow_engine=None, 
              flow_client=None, 
              job_client=None,
              job_dir_factory=None,
              job_runner=None,
              tick_ctx=None,
              flow_runner=None,
              **kwargs
             ):
        self.flow_generator_classes = flow_generator_classes or \
                self.generate_flow_generator_classes()
        self.flow_engine = flow_engine or self.generate_flow_engine(
            flow_generator_classes=self.flow_generator_classes)
        self.flow_client = flow_client or self.generate_flow_client()
        self.job_client = job_client or self.generate_job_client()
        self.job_dir_factory = job_dir_factory or \
                self.generate_job_dir_factory()
        self.job_runner = job_runner or self.generate_job_runner()
        self.tick_ctx = self.decorate_tick_ctx(tick_ctx=tick_ctx)
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

    def generate_flow_client(self):
        return FlowClient(base_url=self.flow_server_url,
                          request_client=self.request_client)

    def generate_job_client(self):
        return JobClient(base_url=self.job_server_url,
                         request_client=self.request_client)

    def generate_job_dir_factory(self):
        return JobDirFactory()

    def generate_job_runner(self): 
        return OdysseyPushJobRunner(job_client=self.job_client,
                                    job_dir_factory=self.job_dir_factory,
                                    odyssey_user=self.odyssey_user,
                                    odyssey_host=self.odyssey_host)

    def decorate_tick_ctx(self, tick_ctx=None):
        tick_ctx = tick_ctx or {}
        decorated_tick_ctx = {
            'create_job': self.job_client.create_job,
            'get_job': self.job_client.fetch_job_by_uuid,
            'create_flow': self.flow_client.create_flow,
            'get_flow': self.flow_client.fetch_flow_by_uuid,
            **tick_ctx
        }
        return decorated_tick_ctx

    def generate_flow_runner(self):
        return FlowRunner(flow_client=self.flow_client,
                          flow_engine=self.flow_engine,
                          tick_ctx=self.tick_ctx)

    def create_flow_record(self, *args, flow_record=None, **kwargs):
        return self.flow_client.create_flow(flow=flow_record)

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
