from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'Run a workflow'

    def add_arguments(self, parser):
        pass
        #parser.add_argument('poll_id', nargs='+', type=int)

    def handle(self, *args, **options):
        self.stdout.write(self.style.SUCCESS('Done.'))
