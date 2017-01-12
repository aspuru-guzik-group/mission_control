import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings
from unittest.mock import MagicMock

from missions.models import Flow as FlowModel
from jobs.models import Job
from flow_runners.base_flow_runner import BaseFlowRunner
from flow_engines.flow_engine import FlowEngine
from flow_client.flow_client import MissionControlFlowClient
from a2g2.flow_generators.reaxys import ReaxysFlowGenerator
from job_runners.base_job_runner import BaseJobRunner
from job_client.job_client import MissionControlJobClient


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
        self.reaxys_flow_spec = self.generate_reaxys_flow_spec()
        self.reaxys_flow_uuid = self.generate_reaxys_flow_model().uuid
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
        return MissionControlJobClient(
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

        def create_flow(flow_kwargs=None):
            flow_spec = flow_kwargs.get('flow_spec', {})
            flow = self.flow_engine.generate_flow(flow_spec=flow_spec)
            serialized_flow = self.flow_engine.serialize_flow(flow)
            created_flow = self.flow_client.create_flow(
                flow={'serialization': json.dumps(serialized_flow)})
            return created_flow['uuid']

        def get_flow(flow_uuid=None):
            return self.flow_client.fetch_flow_by_uuid(uuid=flow_uuid)

        def get_job(job_uuid=None):
            return self.job_client.fetch_job_by_uuid(uuid=job_uuid)

        return {
            'create_flow': create_flow,
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
            job_client=self.job_client,
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
        flow = ReaxysFlowGenerator.generate_flow(
            flow_spec=self.reaxys_flow_spec)
        serialization = self.flow_engine.serialize_flow(flow=flow)
        flow_model = FlowModel.objects.create(
            serialization=json.dumps(serialization))
        return flow_model

    def test_flow(self):
        # flow_runner.tick = 1
        self.flow_runner.tick()
        self.capture_state()
        reaxys_flow = self.states[-1]['flows'][self.reaxys_flow_uuid]
        confgen_task = reaxys_flow.tasks['confgen']
        confgen_flow_uuid = confgen_task.data['flow_uuid']

        self.assertEqual(
            self.states[-1]['flow_models'][confgen_flow_uuid]['status'],
            'PENDING')

        # flow_runner.tick = 2
        self.flow_runner.tick()
        self.capture_state()
        self.assertEqual(
            self.states[-1]['flow_models'][confgen_flow_uuid]['status'],
            'RUNNING')

        # Fake completion of confgen flow.
        FlowModel.objects.filter(uuid=confgen_flow_uuid)\
                .update(status='COMPLETED')

        # flow_runner.tick = 3
        self.flow_runner.tick()
        self.capture_state()
        self.assertEqual(
            self.states[-1]['flow_models'][confgen_flow_uuid]['status'],
            'COMPLETED')
        self.assertEqual(
            self.states[-1]['flow_models'][self.reaxys_flow_uuid]['status'],
            'COMPLETED')

    def capture_state(self):
        if len(self.states) > 1: prev_state = self.states[-1]
        else: prev_state = None
        self.states.append(self.get_state(prev_state=prev_state))

    def get_state(self, prev_state=None):
        state = {
            'flow_models': {
                flow_model.uuid: flow_model.__dict__
                for flow_model in FlowModel.objects.all()
            },
            'jobs': {
                job.uuid: job.__dict__
                for job in Job.objects.all()
            }
        }
        state['flows'] = {
            flow_model['uuid']: self.flow_engine.deserialize_flow(
                serialized_flow=json.loads(flow_model['serialization']))
            for flow_model in state['flow_models'].values()
        }
        return state
