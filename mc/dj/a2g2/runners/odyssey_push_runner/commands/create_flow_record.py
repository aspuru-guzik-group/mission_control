import json

from .base import BaseCommand


class Command(BaseCommand):
    def add_arguments(self, parser=None):
        self.add_flow_spec_arguments(parser=parser)
        self.add_job_and_flow_server_url_arguments(parser)

    def add_flow_spec_arguments(self, parser=None):
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            '--flow_spec_json',
            help="flow_spec as json string.",
        )
        group.add_argument(
            '--flow_spec_file',
            help="path to flow_spec json file.",
        )

    def handle(self, *args, **options):
        self.runner = self.generate_runner_from_options(options=options)
        flow_record = self.assemble_flow_record_from_options(options=options)
        flow = self.runner.create_flow_record(flow_record=flow_record)
        result = {'status': 'SUCCESS', 'flow': flow}
        print(json.dumps(result), file=self.stdout)

    def assemble_flow_record_from_options(self, options):
        if options.get('flow_spec_json', None):
            flow_spec_json = options['flow_spec_json']
        elif options.get('flow_spec_file', None):
            with open(options['flow_spec_file']) as f: flow_spec_json = f.read()
        else: flow_spec_json = None
        flow_record = {'spec': flow_spec_json}
        return flow_record

if __name__ == '__main__':
    Command.run()
