from django.core.management.base import BaseCommand

from ...utils import start_workflow


class Command(BaseCommand):
    help = 'Run a workflow'

    def add_arguments(self, parser):
        parser.add_argument('--spec_key', 
                            help='key for workflow spec to use',
                            required=True, dest='spec_key')

    def handle(self, *args, **options):
        spec_key = options['spec_key']
        start_workflow(spec_key=spec_key)
        self.stdout.write(self.style.SUCCESS('Done.'))
