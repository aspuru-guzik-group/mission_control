from django.core.management.base import BaseCommand

from ...utils import create_spec


class Command(BaseCommand):
    help = 'Create a workflow spec from a serialized format.'

    def add_arguments(self, parser):
        parser.add_argument('--serialization', required=True,
                            help='path to serialized spec')
        parser.add_argument('--serialization_format', default='xml',
                            help='serialization format')

    def handle(self, *args, **options):
        spec = create_spec(**options)
        output_msg = "Created WorkflowSpec '%s'." % spec
        self.stdout.write(self.style.SUCCESS(output_msg))
