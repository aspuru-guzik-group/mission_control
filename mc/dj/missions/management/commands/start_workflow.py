from django.core.management.base import BaseCommand

from ...utils import start_workflow


class Command(BaseCommand):
    help = 'Run a workflow'

    def handle(self, *args, **options):
        start_workflow()
        self.stdout.write(self.style.SUCCESS('Done.'))
