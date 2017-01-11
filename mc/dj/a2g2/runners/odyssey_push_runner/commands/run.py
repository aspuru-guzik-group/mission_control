from .base import BaseCommand
from .. import odyssey_push_runner

class Command(BaseCommand):
    def handle(self, *args, **options):
        self.runner = odyssey_push_runner.OdysseyPushRunner()
        self.runner.run()
