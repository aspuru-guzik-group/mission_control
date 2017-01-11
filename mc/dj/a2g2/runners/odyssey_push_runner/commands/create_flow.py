import json

from .base import BaseCommand
from .. import odyssey_push_runner

class Command(BaseCommand):
    def add_arguments(self, parser=None):
        parser.add_argument(
            '--flow_spec',
            help="flow_spec as json string.",
            type=lambda i: json.loads(i)
        )

    def handle(self, *args, **options):
        self.runner = odyssey_push_runner.OdysseyPushRunner()
        flow = {'spec': options.get('flow_spec')}
        flow_uuid = self.runner.create_flow(flow=flow)
        print("Created flow with uuid '%s'" % flow_uuid, file=self.stdout)
