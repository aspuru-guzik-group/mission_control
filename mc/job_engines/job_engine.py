import json
import importlib
import logging
import os
import random
import tempfile

import dill

from . import constants as _job_engine_constants
from . import utils as _job_engine_utils


class JobEngine(object):
    class JobModuleImportError(Exception): pass

    SUBMISSION_META_NAME = _job_engine_constants.SUBMISSION_META_NAME

    def __init__(self, job_module_loader=None, logger=None):
        self.logger = logger or logging
        self.job_module_loader = job_module_loader or \
                self.get_default_job_module_loader()

    def get_default_job_module_loader(self):
        return DefaultJobModuleLoader()

    def execute_command(self, *args, command=None, **kwargs):
        handler = getattr(self, command)
        return handler(*args, **kwargs)

    def build_job_submission(self, *args, job=None, cfg=None, output_dir=None, 
                             **kwargs):
        output_dir = output_dir or tempfile.mkdtemp()
        os.makedirs(output_dir, exist_ok=True)
        job_module = self.get_job_module(job=job, cfg=cfg)
        submission_meta_from_module = job_module.build_job_submission(
            *args, job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
        submission_meta = {'job': job, 'cfg': cfg, 'dir': output_dir,
                           **submission_meta_from_module}
        self.write_submission_meta(submission_meta=submission_meta,
                                   dir_=output_dir)
        return submission_meta

    def write_submission_meta(self, submission_meta=None, dir_=None):
        submission_meta_path = os.path.join(dir_, self.SUBMISSION_META_NAME)
        with open(submission_meta_path, 'w') as f: json.dump(submission_meta, f)

    def get_job_module(self, job=None, cfg=None):
        try:
            return self.job_module_loader.load_job_module(job=job, cfg=cfg)
        except Exception as exc:
            msg = "Could not load module for job"
            raise self.JobModuleImportError(msg) from exc

    def run_job_submission(self, *args, submission_dir=None, **kwargs):
        submission_meta = self.read_submission_meta(dir_=submission_dir)
        job_module = self.get_job_module(job=submission_meta['job'],
                                         cfg=submission_meta['cfg'])
        return job_module.run_job_submission(
            *args, submission_meta=submission_meta, **kwargs)

    def read_submission_meta(self, dir_=None):
        return _job_engine_utils.read_submission_meta(
            submission_dir=dir_, submission_meta_name=self.SUBMISSION_META_NAME)

class DefaultJobModuleLoader(object):
    def __init__(self):
        self.overrides = {}

    def load_job_module(self, job=None, cfg=None):
        job_module_name = self.get_job_module_name(job=job, cfg=cfg)
        override = self.overrides.get(job_module_name)
        if override:
            module = self.load_module_per_override(
                job=job, cfg=cfg, override=override)
        else:
            module = self.load_module_from_module_path(
                module_path=job_module_name)
        return module

    def load_module_per_override(self, job=None, cfg=None, override=None):
        if isinstance(override, str):
            override = {'type': 'py_module_path', 'params': {'path': override}}
        if override['type'] == 'py_file':
            module = self.load_module_from_file_path(
                file_path=override['params']['path'])
        elif override['type'] == 'py_module':
            module = self.load_module_from_module_path(
                module_path=override['params']['path'])
        elif override['type'] == 'py_obj':
            module = override['params']['obj']
        elif override['type'] == 'dill_file':
            module = self.load_module_from_dill_path(
                dill_path=override['params']['path'])
        return module

    def load_module_from_file_path(self, file_path=None, module_name=None):
        module_name = module_name or 'random_%s' % random.randint(1, int(1e4))
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def load_module_from_module_path(self, module_path=None):
        return importlib.import_module(module_path)

    def load_module_from_dill_path(self, dill_path=None):
        with open(dill_path, 'rb') as f: return dill.load(f)

    def get_job_module_name(self, job=None, cfg=None):
        module_name = job['job_spec']['job_type']
        return module_name
