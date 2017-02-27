import argparse
import json
import os
import sys

class BaseCommand(object):
    help = 'base_command'

    def __init__(self, streams=None):
        self.setup_streams(streams=streams)

    def execute(self, argv=None):
        argv = argv or sys.argv[1:]
        parser = argparse.ArgumentParser()
        self.add_arguments(parser)
        parsed_args = parser.parse_args(argv)
        self.handle(**vars(parsed_args))

    def setup_streams(self, streams=None):
        streams = streams or {}
        for stream_id in ['stdout', 'stderr', 'stdin']:
            if stream_id in streams:
                stream = streams[stream_id]
            else:
                stream = getattr(sys, stream_id)
            setattr(self, stream_id, stream)

    def add_arguments(self, parser):
        def json_file_type(file_path): return json.load(open(file_path))
        parser.add_argument('--job', type=json_file_type)
        parser.add_argument('--cfg', type=json_file_type, default={})
        parser.add_argument('--ctx_dir', type=str)
        parser.add_argument('--output_dir', type=str)

    def get_cfg_value(self, cfg=None, key=None, default=None):
        if key in os.environ: cfg_value = os.environ[key]
        else: cfg_value = (cfg or {}).get(key, default)
        return cfg_value

    def handle(self, *args, **kwargs):
        self.execute_job(*args, **kwargs)

    def execute_job(self, *args, job=None, cfg=None, ctx_dir=None, 
                    output_dir=None, **kwargs): pass
