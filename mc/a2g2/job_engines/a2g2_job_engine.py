import argparse
import importlib
import json
import sys

DEFAULT_JOB_MODULE_PKGS = ['mc.a2g2.job_modules']
DEFAULT_JOB_TYPE_PREFIX = 'a2g2.jobs.'

class A2G2JobEngine(object):
    """Facade for job lifecycle commands and job module dispatch."""

    def execute_command(self, *args, command=None, **kwargs):
        handler = getattr(self, command)
        handler(*args, **kwargs)

    def build_job_submission(self, *args, job=None, cfg=None, job_dir=None, 
                             **kwargs):
        job_module = self.get_job_module(job=job, cfg=cfg)
        return job_module.build_job_submission(
            *args, job=job, cfg=cfg, job_dir=job_dir, **kwargs)

    def run_job_submission(self, *args, job=None, cfg=None, submission=None,
                           **kwargs):
        job_module = self.get_job_module(job=job, cfg=cfg)
        return job_module.run_job_submission(
            *args, job=job, cfg=cfg, submision=submission, **kwargs)

    def get_job_module(self, job=None, cfg=None):
        job_module_pkgs = DEFAULT_JOB_MODULE_PKGS
        job_module_name = self.get_job_module_name(job=job, cfg=cfg)
        try:
            for job_module_pkg in job_module_pkgs:
                full_module_name = '.'.join([job_module_pkg, job_module_name])
                try: return importlib.import_module(full_module_name)
                except ImportError: pass
            no_module_found_msg = ("No job_module found for job_module"
                                   " '{job_module_name}'").format(
                                       job_module_name=job_module_name)
            raise Exception(no_module_found_msg)
        except Exception as error:
            raise Exception("Could not load module for job: %s" % error)

    def get_job_module_name(self, job=None, cfg=None):
        module_name = self.get_module_name_from_cfg(job=job, cfg=cfg) or \
                self.get_module_name_from_type(job=job, cfg=cfg)
        return module_name

    def get_module_name_from_cfg(self, job=None, cfg=None):
        engine_cfg = cfg.get('job_engine', {})
        modules_for_types = engine_cfg.get('modules_for_types', {})
        return modules_for_types.get(job['job_spec']['job_type'], None)

    def get_module_name_from_type(self, job=None, cfg=None):
        engine_cfg = cfg.get('job_engine', {})
        job_type_prefix = engine_cfg.get(
            'job_type_prefix', DEFAULT_JOB_TYPE_PREFIX)
        return job['job_spec']['job_type'].replace(job_type_prefix, '')

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
        def json_file_type(file_path): return json.load(open(file_path))
        parser.add_argument('--job', type=json_file_type)
        parser.add_argument('--cfg', type=json_file_type, default={})
        subparsers = parser.add_subparsers(dest='command',
                                           title='commands',
                                           description='valid commands')
        subparsers.add_parser('build_job_submission')
        run_job_submission_parser = subparsers.add_parser('run_job_submission')
        run_job_submission_parser.add_argument('--submission',
                                               type=json_file_type)

    def handle(self, *args, command=None, job=None, cfg=None, **kwargs):
        self.execute_command(job_engine=self.generate_job_engine(cfg=cfg),
                             command=command, job=job, cfg=cfg, **kwargs)

    def generate_job_engine(self, cfg=None):
        cfg = cfg or {}
        engine_cfg = cfg.get('job_engine', {})
        # Class spec allows for using a mock engine in tests.
        class_spec = engine_cfg.get('engine_class_spec', None)
        if class_spec:
            JobEngineClass = self.get_class_from_class_spec(class_spec)
        else:
            JobEngineClass = A2G2JobEngine
        job_engine = JobEngineClass()
        return job_engine

    def get_class_from_class_spec(self, class_spec=None):
        class_spec_parts = class_spec.split('.')
        module_name = '.'.join(class_spec_parts[:-1])
        class_name = class_spec_parts[-1]
        job_module = importlib.import_module(module_name)
        Clazz = getattr(job_module, class_name)
        return Clazz

    def execute_command(self, job_engine=None, command=None, job=None, cfg=None,
                        **kwargs):
        return job_engine.execute_command(command=command, job=job, cfg=cfg,
                                          **kwargs)

if __name__ == '__main__':
    JobEngineCommand().execute(argv=sys.argv[1:])
