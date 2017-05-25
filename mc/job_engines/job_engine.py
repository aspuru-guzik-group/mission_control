import json
import importlib
import logging
import os
import tempfile


class JobEngine(object):
    class JobModuleImportError(Exception): pass

    SUBMISSION_META_NAME = 'MC_JOB_ENGINE__SUBMISSION_META.json'

    def __init__(self, job_module_loader=None, logger=None):
        self.logger = logger or logging
        self.job_module_loader = job_module_loader or \
                self.get_default_job_module_loader()

    def get_default_job_module_loader(self):
        return DefaultJobModuleLoader()

    def execute_command(self, *args, command=None, **kwargs):
        handler = getattr(self, command)
        return handler(*args, **kwargs)

    def build_submission(self, *args, job=None, cfg=None, output_dir=None, 
                         **kwargs):
        output_dir = output_dir or tempfile.mkdtemp()
        os.makedirs(output_dir, exist_ok=True)
        job_module = self.get_job_module(job=job, cfg=cfg)
        submission_meta_from_module = job_module.build_submission(
            *args, job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
        submission_meta = {'job': job, 'cfg': cfg, 'dir': output_dir,
                           **submission_meta_from_module}
        self.write_submission_meta(submission_meta=submission_meta,
                                   _dir=output_dir)
        return submission_meta

    def write_submission_meta(self, submission_meta=None, _dir=None):
        submission_meta_path = os.path.join(_dir, self.SUBMISSION_META_NAME)
        with open(submission_meta_path, 'w') as f: json.dump(submission_meta, f)

    def get_job_module(self, job=None, cfg=None):
        try:
            return self.job_module_loader.load_job_module(job=job, cfg=cfg)
        except Exception as exc:
            msg = "Could not load module for job"
            raise self.JobModuleImportError(msg) from exc

    def run_submission(self, *args, submission_dir=None, **kwargs):
        submission_meta = self.read_submission_meta(_dir=submission_dir)
        job_module = self.get_job_module(job=submission_meta['job'],
                                         cfg=submission_meta['cfg'])
        return job_module.run_submission(
            *args, submission_meta=submission_meta, **kwargs)

    def read_submission_meta(self, _dir=None):
        submission_meta_path = os.path.join(_dir, self.SUBMISSION_META_NAME)
        with open(submission_meta_path) as f: return json.load(f)

class DefaultJobModuleLoader(object):
    def load_job_module(self, job=None, cfg=None):
        job_module_name = self.get_job_module_name(job=job, cfg=cfg)
        return importlib.import_module(job_module_name)

    def get_job_module_name(self, job=None, cfg=None):
        module_name = job['job_spec']['job_type']
        return module_name
