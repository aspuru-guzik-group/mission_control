import json

from django.core.management.base import BaseCommand

from ...models import Job

class Command(BaseCommand):
    help = 'Create a job.'

    def add_arguments(self, parser):
        parser.add_argument('--spec', 
                            help='path to spec to use for job')

    def handle(self, *args, **options):
        if options['spec']:
            spec = self.load_spec(spec_path=options['spec'])
        else:
            spec = {}
        job = Job.objects.create(spec=spec)
        self.stdout.write("Created job '%s'." % job)

    def load_spec(self, spec_path=None):
        return json.load(open(spec_path))
