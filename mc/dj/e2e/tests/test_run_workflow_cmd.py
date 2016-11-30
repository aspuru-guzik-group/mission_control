from django.conf.urls import url, include
from django.core.management import call_command
from django.test import TestCase, override_settings
import json

from jobs.models import Job, JobStatuses
from .stub_runner import generate_stub_runner
from missions.models import Workflow

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

    def test_run_workflow(self):
        call_command('run_workflow')

        workflows = Workflow.objects.all()
        self.assertEqual(len(workflows), 1)

        # First workflow tick: GenerateNames.
        jobs = Job.objects.order_by('modified')
        self.assertEqual(jobs[0].type, 'GenerateNames')
        self.job_runner.run(ntimes=3)

        # Second workflow tick: IngestNames.
        jobs = Job.objects.order_by('modified')
        self.assertEqual(len(jobs), 2)
        self.assertEqual(jobs[1].type, 'IngestNames')
        self.assertEqual(jobs[1].status, JobStatuses.Pending.name)
        self.job_runner.run(ntimes=3)
        self.assertEqual(jobs[1].status, JobStatuses.Completed.name)

        # Finalized.
        workflows = Workflow.objects.all()
        self.assertEqual(workflows[0].finished, True)
