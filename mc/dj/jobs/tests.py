from django.test import TestCase
from django.test.utils import override_settings
from django.urls import reverse
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
        self.assertTrue(job.uuid is not None)
        self.assertTrue(job.created is not None)
        self.assertTrue(job.modified is not None)

@override_settings(ROOT_URLCONF='jobs.urls')
class JobsAPITestCase(APITestCase):
    def setUp(self):
        self.jobs = [Job.objects.create(name="job_%s" % i)
                     for i in range(1)]

    def test_list_jobs(self):
        response = self.client.get('/jobs/')
        expected_data = [JobSerializer(job).data for job in self.jobs]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))
