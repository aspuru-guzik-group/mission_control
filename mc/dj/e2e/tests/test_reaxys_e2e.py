import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings
from unittest.mock import MagicMock

from missions.models import Flow as FlowModel
from jobs.models import Job
from flow_runners.base_flow_runner import BaseFlowRunner
from flow_engines.flow_engine import FlowEngine
from flow_client.flow_client import MissionControlFlowClient
from a2g2_flows.flow_generators.reaxys import ReaxysFlowGenerator
from job_runners.base_job_runner import BaseJobRunner
from job_spec_client.job_spec_client import MissionControlJobSpecClient


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('missions.urls')),
    url(r'^%s' % BASE_PATH, include('jobs.urls')),
]

@override_settings(ROOT_URLCONF=__name__)
class ReaxysFlowE2ETestCase(TestCase):
    def setUp(self):
        self.configure_request_client()
        self.flow_client = self.generate_flow_client()
        self.job_client = self.generate_job_client()
        self.flow_engine = self.generate_flow_engine()
        self.tick_ctx = self.generate_tick_ctx()
        self.flow_runner = self.generate_flow_runner()
        self.job_runner = self.generate_job_runner()
        self.flow_spec = self.generate_reaxys_flow_spec()
        self.flow_model = self.generate_reaxys_flow_model()
        self.states = []

    def configure_request_client(self):
        orig_patch = self.client.patch
        def json_patch(path, data=None, **kwargs):
            return orig_patch(path, json.dumps(data),
                              content_type='application/json', **kwargs)
        self.client.patch = json_patch

    def generate_flow_client(self):
        return MissionControlFlowClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )

    def generate_job_client(self):
        return MissionControlJobSpecClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )

    def generate_flow_engine(self):
        flow_engine = FlowEngine()
        for flow_generator_class in [ReaxysFlowGenerator]:
            flow_engine.register_flow_generator_class(
                flow_generator_class=flow_generator_class)
        return flow_engine

    def generate_tick_ctx(self):
        def get_flow(flow_id=None):
            raise NotImplementedError

        def get_job(job_id=None):
            raise NotImplementedError

        return {
            'create_flow': self.flow_client.create_flow,
            'get_flow': get_flow,
            'create_job': self.job_client.create_job,
            'get_job': get_job,
        }

    def generate_flow_runner(self):
        return BaseFlowRunner(flow_client=self.flow_client,
                              flow_engine=self.flow_engine,
                              tick_ctx=self.tick_ctx)

    def generate_job_runner(self):
        execution_client = MagicMock()
        job_dir_factory = MagicMock()
        transfer_client = MagicMock()
        runner = BaseJobRunner(
            job_spec_client=self.job_client,
            execution_client=execution_client,
            job_dir_factory=job_dir_factory,
            transfer_client=transfer_client,
        )
        return runner

    def generate_reaxys_flow_spec(self):
        return {
            'flow_type': ReaxysFlowGenerator.flow_type,
            'smiles': 'some smiles',
            'confgen_params': 'some confgen params',
        }

    def generate_reaxys_flow_model(self):
        flow = ReaxysFlowGenerator.generate_flow(flow_spec=self.flow_spec)
        serialization = self.flow_engine.serialize_flow(flow=flow)
        flow_model = FlowModel.objects.create(
            serialization=json.dumps(serialization))
        return flow_model

    def test_flow(self):
        states = []
        states.append(self.get_state(prev_state=None))
        self.tick_and_capture_state()
        self.fail()

    def tick_and_capture_state(self):
        self.flow_runner.tick()
        self.states.append(self.get_state(prev_state=self.states[-1]))

    def get_state(self, prev_state=None):
        state = {
            'flows': {
                flow_model.uuid: flow_model.__dict__
                for flow_model in FlowModel.objects.all()
            },
            'jobs': {
                job.uuid: job.__dict__
                for job in Job.objects.all()
            }
        }
        return state
