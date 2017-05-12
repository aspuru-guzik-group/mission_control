import json
import importlib
import os


class JobEngine(object):
    class JobModuleImportError(Exception): pass

    SUBMISSION_META_NAME = 'MC_JOB_ENGINE__SUBMISSION_META.json'

    def execute_command(self, *args, command=None, **kwargs):
        handler = getattr(self, command)
        return handler(*args, **kwargs)

    def build_submission(self, *args, job=None, cfg=None, output_dir=None, 
                             **kwargs):
        job_module = self.get_job_module(job=job, cfg=cfg)
        submission_meta_from_module = job_module.build_submission(
            *args, job=job, cfg=cfg, output_dir=output_dir, **kwargs)
        submission_meta = {'job': job, 'cfg': cfg,
                           **submission_meta_from_module}
        self.write_submission_meta(submission_meta=submission_meta,
                                   _dir=output_dir)
        return submission_meta

    def write_submission_meta(self, submission_meta=None, _dir=None):
        submission_meta_path = os.path.join(_dir, self.SUBMISSION_META_NAME)
        with open(submission_meta_path, 'w') as f: json.dump(submission_meta, f)

    def get_job_module(self, job=None, cfg=None):
        try:
            job_module_name = self.get_job_module_name(job=job, cfg=cfg)
            return importlib.import_module(job_module_name)
        except Exception as exc:
            raise self.JobModuleImportError(
                "Could not load module for job") from exc

    def get_job_module_name(self, job=None, cfg=None):
        module_name = job['job_spec']['job_type']
        return module_name

    def run_submission(self, *args, submission_dir=None, **kwargs):
        submission_meta = self.read_submission_meta(_dir=submission_dir)
        job_module = self.get_job_module(job=submission_meta['job'],
                                         cfg=submission_meta['cfg'])
        return job_module.run_submission(
            *args, submission_meta=submission_meta, **kwargs)

    def read_submission_meta(self, _dir=None):
        submission_meta_path = os.path.join(_dir, self.SUBMISSION_META_NAME)
        with open(submission_meta_path) as f: return json.load(f)
