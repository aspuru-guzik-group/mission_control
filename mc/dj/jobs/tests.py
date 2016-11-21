import json
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from .models import Job
from .serializers import JobSerializer

class JobTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'name': 'test_name',
        }
        job = Job.objects.create(**kwargs)
        self.assertEqual(job.name, kwargs['name'])
        self.assertEqual(job.status, Job.STATUSES.PENDING.name)
        self.assertTrue(job.uuid is not None)
        self.assertTrue(job.created is not None)
        self.assertTrue(job.modified is not None)

@override_settings(ROOT_URLCONF='jobs.urls')
class ListJobsTestCase(APITestCase):
    def setUp(self):
        self.jobs = [Job.objects.create(name="job_%s" % i)
                     for i in range(1)]

    def test_list_jobs(self):
        response = self.client.get('/jobs/')
        expected_data = [JobSerializer(job).data for job in self.jobs]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

@override_settings(ROOT_URLCONF='jobs.urls')
class PatchJobTestCase(APITestCase):
    def setUp(self):
        self.jobs = [Job.objects.create(name="job_%s" % i)
                     for i in range(1)]

    def test_patch_job(self):
        job_to_patch = self.jobs[0]
        new_values = {'name': 'new_name'}
        response = self.client.patch('/jobs/%s/' % job_to_patch.uuid,
                                     new_values)
        self.assertEqual(response.status_code, 200)
        patched_job = Job.objects.get(uuid=job_to_patch.uuid)
        patched_job_attrs = {attr: getattr(patched_job, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_job_attrs, new_values)

    def test_patch_with_bad_status(self):
        response = self.client.patch('/jobs/%s/' % self.jobs[0].uuid,
                                     {'status': 'BADDO!'})
        self.assertEqual(response.status_code, 400)

@override_settings(ROOT_URLCONF='jobs.urls')
class ClaimJobTestCase(TestCase):
    def setUp(self):
        self.unclaimed_jobs = [Job.objects.create(name="job_%s" % i)
                               for i in range(1)]
        self.claimed_jobs = [Job.objects.create(
            name="job_%s" % i, status=Job.STATUSES.CLAIMED.name)
            for i in range(1)]
        self.all_jobs = self.unclaimed_jobs + self.claimed_jobs
        self.jobs_to_claim = [j for j in self.unclaimed_jobs[:-1]] + [
            j for j in self.claimed_jobs[:-1]]
        self._claimJobs()

    def _claimJobs(self):
        csv_uuids = ','.join([str(job.uuid) for job in self.jobs_to_claim])
        self.response = self.client.post('/claim_jobs/', {'uuids': csv_uuids})

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {str(job.uuid): (job in self.unclaimed_jobs)
                         for job in self.jobs_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        statuses = {job.uuid: job.status for job in Job.objects.all()}
        expected_statuses = {}
        for job in self.all_jobs:
            if job in set(self.jobs_to_claim + self.claimed_jobs):
                expected_status = Job.STATUSES.CLAIMED.name
            else:
                expected_status = Job.STATUSES.PENDING.name
            expected_statuses[job.uuid] = expected_status
        self.assertEqual(statuses, expected_statuses)

