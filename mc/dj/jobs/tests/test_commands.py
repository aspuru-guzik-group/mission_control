import json
import os
import tempfile
import io

from django.core.management import call_command
from django.test import TestCase

from ..models import Job


class CreateJobTestCase(TestCase):
    def setUp(self):
        self.setup_spec()
        self.cmd_output = self.call_create_job_command(spec_path=self.spec_path)
        self.created_job = Job.objects.all()[0]

    def setup_spec(self):
        self.spec = {'some': 'spec'}
        tmp_dir = tempfile.mkdtemp()
        self.spec_path = os.path.join(tmp_dir, 'spec.json')
        json.dump(self.spec, open(self.spec_path, 'w'))

    def test_creates_job_with_expected_attrs(self):
        expected_attrs = {'spec': self.spec}
        actual_attrs = {attr: getattr(self.created_job, attr)
                        for attr in expected_attrs}
        self.assertEqual(actual_attrs, expected_attrs)

    def test_produces_expected_output(self):
        self.assertEqual(self.cmd_output.strip(),
                         "Created job '%s'." % self.created_job)

    def call_create_job_command(self, spec_path=None):
        f = io.StringIO()
        args = ['--spec', spec_path]
        call_command('create_job', *args, stdout=f)
        return f.getvalue()
