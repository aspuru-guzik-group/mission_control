import argparse
import json
import sys

import requests as request_client

from .. import odyssey_push_runner

class BaseCommand(object):
    @classmethod
    def run(cls, args=sys.argv, streams=None):
        command = cls()
        command.set_streams(**(streams or {}))
        parser = argparse.ArgumentParser()
        command.add_arguments(parser)
        command.add_params_file_argument(parser)
        options = vars(parser.parse_args(args=args))
        params_file_path = options.get('params_file', None)
        if params_file_path:
            params_from_file = json.load(open(params_file_path))
            options.update(params_from_file)
            del options['params_file']
        command.handle(**options)

    def __init__(self, *args, **kwargs):
        self.set_streams(**kwargs)

    def set_streams(self, **kwargs):
        for stream_id in ['stdin', 'stdout', 'stderr']:
            stream = kwargs.get(stream_id, getattr(sys, stream_id))
            setattr(self, stream_id, stream)

    def add_arguments(self, parser): pass

    def add_params_file_argument(self, parser):
        parser.add_argument(
            '--params_file',
            help=("JSON file containing params. Any parameters in this file"
                  " will override parameters provided on the command line."),
        )

    def add_job_and_flow_server_url_arguments(self, parser):
        self.add_job_server_url_argument(parser)
        self.add_flow_server_url_argument(parser)

    def add_job_server_url_argument(self, parser):
        parser.add_argument('--job_server_url', help="URL of job server")

    def add_flow_server_url_argument(self, parser):
        parser.add_argument('--flow_server_url', help="URL of flow server")

    def handle(self, *args, **kwargs): pass

    def generate_runner_from_options(self, options=None):
        return odyssey_push_runner.OdysseyPushRunner(
            request_client=request_client, **options)
