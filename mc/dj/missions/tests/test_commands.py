import json
import os
import tempfile
import io

from django.core.management import call_command
from django.test import TestCase

from ..models import Job


class CreateJobTestCase(TestCase):
    def setUp(self):
        self.setup_job_spec()
        self.cmd_output = self.call_create_job_command(job_spec_path=self.job_spec_path)
        self.created_job = Job.objects.all()[0]

    def setup_job_spec(self):
        self.job_spec = {'some': 'job_spec'}
        tmp_dir = tempfile.mkdtemp()
        self.job_spec_path = os.path.join(tmp_dir, 'job_spec.json')
        json.dump(self.job_spec, open(self.job_spec_path, 'w'))

    def test_creates_job_with_expected_attrs(self):
        expected_attrs = {'job_spec': self.job_spec}
        actual_attrs = {attr: getattr(self.created_job, attr)
                        for attr in expected_attrs}
        self.assertEqual(actual_attrs, expected_attrs)

    def test_produces_expected_output(self):
        self.assertEqual(self.cmd_output.strip(),
                         "Created job '%s'." % self.created_job)

    def call_create_job_command(self, job_spec_path=None):
        f = io.StringIO()
        args = ['--job_spec', job_spec_path]
        call_command('create_job', *args, stdout=f)
        return f.getvalue()
