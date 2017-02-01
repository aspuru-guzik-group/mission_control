import json

from django.core.management.base import BaseCommand

from ...models import Job

class Command(BaseCommand):
    help = 'Create a job.'

    def add_arguments(self, parser):
        parser.add_argument('--job_spec', 
                            help='path to spec to use for job')

    def handle(self, *args, **options):
        if options['job_spec']:
            job_spec = self.load_job_spec(spec_path=options['job_spec'])
        else:
            job_spec = {}
        job = Job.objects.create(job_spec=job_spec)
        self.stdout.write("Created job '%s'." % job)

    def load_job_spec(self, job_spec_path=None):
        return json.load(open(job_spec_path))
