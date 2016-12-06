import json

from django.conf.urls import url, include
from django.core.management import call_command
from django.test import TestCase, override_settings
import SpiffWorkflow.specs

from jobs.models import Job, JobStatuses
from .stub_runner import generate_stub_runner
from missions.models import Workflow, WorkflowSpec, WorkflowStatuses
from missions import utils as mission_utils

BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('jobs.urls')),
]

@override_settings(ROOT_URLCONF=__name__)
class JobRunner_X_Workflow_TestCase(TestCase):
    def setUp(self):
        self.job_runner = generate_stub_runner(base_url='/%s' % BASE_PATH,
                                               request_client=self.client)
        self.job_runner.tick_interval = .1
        orig_patch = self.client.patch
        def json_patch(path, data=None, **kwargs):
            return orig_patch(path, json.dumps(data),
                              content_type='application/json', **kwargs)
        self.client.patch = json_patch
        self.spec = self.generate_workflow_spec()

    def generate_workflow_spec(self):
        spec = WorkflowSpec.objects.create(
            key='my_spec_key',
            serialization=self.generate_serialized_spiff_spec()
        )
        return spec

    def generate_serialized_spiff_spec(self):
        spiff_spec = self.generate_spiff_spec()
        workflow_engine = mission_utils.get_workflow_engine()
        return workflow_engine.serialize_spiff_spec(spiff_spec=spiff_spec)

    def generate_spiff_spec(self):
        spec = SpiffWorkflow.specs.WorkflowSpec()
        tasks = {
            task_name: SpiffWorkflow.specs.Simple(spec, task_name, data={
                'job_spec': {'type': task_name}
            }) for task_name in ['GenerateNames', 'IngestNames']
        }
        spec.start.connect(tasks['GenerateNames'])
        tasks['GenerateNames'].connect(tasks['IngestNames'])
        return spec

    def test_start_workflow(self):
        call_command('start_workflow', *['--spec_key', self.spec.key])

        workflows = Workflow.objects.all()
        self.assertEqual(len(workflows), 1)

        # First workflow tick: GenerateNames.
        jobs = Job.objects.order_by('modified')
        self.assertEqual(jobs[0].spec['type'], 'GenerateNames')
        self.job_runner.run(ntimes=3)

        # Second workflow tick: IngestNames.
        mission_utils.poll_workflow_jobs()
        jobs = Job.objects.order_by('modified')
        self.assertEqual(len(jobs), 2)
        self.assertEqual(jobs[1].spec['type'], 'IngestNames')
        self.assertEqual(jobs[1].status, JobStatuses.Pending.name)
        self.job_runner.run(ntimes=3)
        jobs = Job.objects.order_by('modified')
        self.assertEqual(jobs[1].status, JobStatuses.Completed.name)

        # Final tick.
        mission_utils.poll_workflow_jobs()
        workflows = Workflow.objects.all()
        self.assertEqual(workflows[0].status, WorkflowStatuses.Completed.name)
