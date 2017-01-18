from django.test import TestCase

from ..models import Job, JobStatuses


class JobTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'name': 'test_name',
        }
        job = Job.objects.create(**kwargs)
        self.assertEqual(job.name, kwargs['name'])
        self.assertEqual(job.status, JobStatuses.PENDING.name)
        expected_attrs = ['uuid', 'created', 'modified', 'spec']
        for attr in expected_attrs:
            self.assertTrue(hasattr(job, attr))
