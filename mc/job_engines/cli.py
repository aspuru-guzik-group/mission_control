import argparse
import importlib
import json
import os
import random
import sys

import dill

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
        parser.add_argument('--engine_config', type=os.path.expanduser,
                            help=("Path to config file."
                                  " File type will be determined by extension,"
                                  " or by --engine_config_type argument"))
        parser.add_argument('--engine_config_type')
        subparsers = parser.add_subparsers(dest='command',
                                           title='commands',
                                           description='valid commands')
        self.add_build_submission_subparser(subparsers=subparsers)
        self.add_run_submission_subparser(subparsers=subparsers)

    def _json_file(self, file_path):
        with open(file_path) as f: return json.load(f)

    def add_build_submission_subparser(self, subparsers=None):
        subparser = subparsers.add_parser('build_submission')
        subparser.add_argument('--job', type=self._json_file)
        subparser.add_argument('--cfg', type=self._json_file, default={})
        subparser.add_argument('--output_dir')

    def add_run_submission_subparser(self, subparsers=None):
        subparser = subparsers.add_parser('run_submission')
        subparser.add_argument('--submission_dir')

    def handle(self, *args, command=None, **kwargs):
        self.execute_command(job_engine=self.get_job_engine(**kwargs),
                             command=command, **kwargs)

    def get_job_engine(self, *args, engine_config=None, engine_config_type=None,
                       **kwargs):
        if engine_config:
            engine_config = self.load_from_file(file_path=engine_config,
                                                file_type=engine_config_type)
            job_engine = getattr(engine_config, 'job_engine')
        else:
            job_engine = JobEngine()
        return job_engine

    def load_from_file(self, file_path=None, file_type=None):
        file_type = file_type or os.path.splitext(file_path)[-1]
        if file_type == '.py':
            return self.load_module_from_file_path(file_path=file_path)
        elif file_type == '.dill':
            with open(file_path, 'rb') as f: return dill.load(f)
        else: raise Exception(f"unknown file_type '{file_type}'")

    def load_module_from_file_path(self, file_path=None, module_name=None):
        module_name = module_name or 'random_%s' % random.randint(1, int(1e4))
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def execute_command(self, job_engine=None, command=None, **kwargs):
        return job_engine.execute_command(command=command, **kwargs)

if __name__ == '__main__':
    JobEngineCommand().execute(argv=sys.argv[1:])
