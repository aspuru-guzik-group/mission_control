import argparse
import json
import sys

class BaseCommand(object):
    @classmethod
    def run(cls, args=sys.argv):
        command = cls()
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

    def handle(self, *args, **kwargs): pass
