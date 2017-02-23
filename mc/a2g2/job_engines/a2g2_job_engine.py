import argparse
import importlib
import json
import sys


class A2G2JobEngine(object):
    def execute_job(self, job=None, cfg=None, output_dir=None):
        job_module = self.get_job_module(job=job, cfg=cfg)
        return job_module.execute_job(job=job, cfg=cfg, output_dir=output_dir)

    def get_job_module(self, job=None, cfg=None):
        try:
            module_name = '{job_type}_job_module'.format(
                job_type=job['job_spec']['job_type'])
            job_module = importlib.import_module(
                module_name, package=sys.modules[__name__].__package__)
        except Exception as error:
            raise Exception("Could not load job_module for job: %s" % error)
        return job_module

class ExecuteJobCommand(object):
    help = 'Execute Job'

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
        parser.add_argument('--output_dir', type=str)

    def handle(self, *args, job=None, cfg=None, output_dir=None, **kwargs):
        self.execute_job(
            job=job,
            cfg=cfg,
            output_dir=output_dir,
            job_engine=self.generate_job_engine(cfg=cfg),
        )

    def generate_job_engine(self, cfg=None):
        cfg = cfg or {}
        class_spec = cfg.get('A2G2_JOB_ENGINE_CLASS', None)
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

    def execute_job(self, job=None, cfg=None, output_dir=None,
                    job_engine=None):
        return job_engine.execute_job(job=job, cfg=cfg, output_dir=output_dir)

if __name__ == '__main__':
    cmd = ExecuteJobCommand()
    cmd.execute(argv=sys.argv[1:])
