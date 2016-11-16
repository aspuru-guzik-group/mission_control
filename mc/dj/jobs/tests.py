from django.test import TestCase

from .models import Job

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
