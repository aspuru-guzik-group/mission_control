from .base import BaseCommand
from .. import odyssey_push_runner

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--odyssey_user', help="Odyssey user")
        parser.add_argument('--odyssey_host', help="Odyssey host")

    def handle(self, *args, **options):
        self.runner = odyssey_push_runner.OdysseyPushRunner(**options)
        self.runner.run()
