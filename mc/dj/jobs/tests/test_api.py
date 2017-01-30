import json
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..models import Job, JobStatuses
from ..serializers import JobSerializer

@override_settings(ROOT_URLCONF='jobs.urls')
class ListJobsTestCase(APITestCase):
    def test_list_jobs(self):
        jobs = [Job.objects.create(name="job_%s" % i) for i in range(3)]
        response = self.client.get('/jobs/')
        expected_data = [JobSerializer(job).data for job in jobs]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

    def test_status_filtering(self):
        statuses = [JobStatuses.PENDING.name, JobStatuses.RUNNING.name]
        jobs_by_status = {
            status: [Job.objects.create(name="job_%s" % i, status=status)
                     for i in range(3)]
            for status in statuses
        }
        for status in statuses:
            response = self.client.get('/jobs/', {'status': status})
            expected_data = [JobSerializer(job).data
                             for job in jobs_by_status[status]]
            self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                             sorted(expected_data, key=lambda j: j['uuid']))

    def test_uuid_filtering(self):
        jobs = [Job.objects.create() for i in range(2)]
        jobs_by_uuid = {job.uuid: job for job in jobs}
        for uuid, job in jobs_by_uuid.items():
            response = self.client.get('/jobs/', {'uuid': uuid})
            expected_data = [JobSerializer(job).data]
            self.assertEqual(self.key_by_uuid(response.data),
                             self.key_by_uuid(expected_data))

    def key_by_uuid(self, items=None):
        return {item['uuid']: item for item in items}

@override_settings(ROOT_URLCONF='jobs.urls')
class PatchJobTestCase(APITestCase):
    def setUp(self):
        self.jobs = [Job.objects.create(name="job_%s" % i)
                     for i in range(1)]

    def test_patch_job(self):
        job_to_patch = self.jobs[0]
        new_values = {'name': 'new_name', 'data': json.dumps({'new': 'data'})}
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
            name="job_%s" % i, status=JobStatuses.RUNNING.name)
            for i in range(1)]
        self.all_jobs = self.unclaimed_jobs + self.claimed_jobs
        self.jobs_to_claim = [j for j in self.unclaimed_jobs[:-1]] + [
            j for j in self.claimed_jobs[:-1]]
        self._claimJobs()

    def _claimJobs(self):
        uuids = [job.uuid for job in self.jobs_to_claim]
        self.response = self.client.post('/claim_jobs/',
                                         json.dumps({'uuids': uuids}),
                                         content_type='application/json')

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {str(job.uuid): JobSerializer(job).data
                         for job in self.jobs_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        statuses = {job.uuid: job.status for job in Job.objects.all()}
        expected_statuses = {}
        for job in self.all_jobs:
            if job in set(self.jobs_to_claim + self.claimed_jobs):
                expected_status = JobStatuses.RUNNING.name
            else:
                expected_status = JobStatuses.PENDING.name
            expected_statuses[str(job.uuid)] = expected_status
        self.assertEqual(statuses, expected_statuses)

