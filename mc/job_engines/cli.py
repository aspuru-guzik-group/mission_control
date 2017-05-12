import argparse
import json
import sys

from .job_engine import JobEngine


class JobEngineCommand(object):
    help = 'JobEngine Command'

    def __init__(self, streams=None):
        self.setup_streams(streams=streams)

    def setup_streams(self, streams=None):
        streams = streams or {}
        for stream_id in ['stdout', 'stderr', 'stdin']:
            if stream_id in streams:
                stream = streams[stream_id]
            else:
                stream = getattr(sys, stream_id)
            setattr(self, stream_id, stream)

    def execute(self, argv=None):
        parser = argparse.ArgumentParser()
        self.add_arguments(parser)
        parsed_args = parser.parse_args(argv)
        self.handle(**vars(parsed_args))

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest='command',
                                           title='commands',
                                           description='valid commands')
        self.add_build_submission_subparser(subparsers=subparsers)
        self.add_run_submission_subparser(subparsers=subparsers)

    def _json_file(self, file_path): return json.load(open(file_path))

    def add_build_submission_subparser(self, subparsers=None):
        subparser = subparsers.add_parser('build_submission')
        subparser.add_argument('--job', type=self._json_file)
        subparser.add_argument('--cfg', type=self._json_file, default={})
        subparser.add_argument('--output_dir')

    def add_run_submission_subparser(self, subparsers=None):
        subparser = subparsers.add_parser('run_submission')
        subparser.add_argument('--submission_dir')

    def handle(self, *args, command=None, **kwargs):
        self.execute_command(job_engine=self.get_job_engine(), command=command,
                             **kwargs)

    def execute_command(self, job_engine=None, command=None, **kwargs):
        return job_engine.execute_command(command=command, **kwargs)

    def get_job_engine(self): return JobEngine()

if __name__ == '__main__':
    JobEngineCommand().execute(argv=sys.argv[1:])
