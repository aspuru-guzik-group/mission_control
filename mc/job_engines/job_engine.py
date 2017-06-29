import json
import logging
import os
import tempfile

from . import constants as _job_engine_constants
from . import utils as _job_engine_utils
from .default_job_module_loader import DefaultJobModuleLoader


class JobEngine(object):
    class JobModuleImportError(Exception): pass

    SUBMISSION_META_NAME = _job_engine_constants.SUBMISSION_META_NAME

    def __init__(self, job_module_loader=None, logger=None):
        """
        Args:
            job_module_loader (JobModuleLoader): loader to use for loading
                job modules. Default: DefaultJobModuleLoader.
            logger: (logging.Logger): a logger.
        """
        self.logger = logger or logging
        self.job_module_loader = job_module_loader or \
                self.get_default_job_module_loader()

    def get_default_job_module_loader(self): return DefaultJobModuleLoader()

    def execute_command(self, command=None, **kwargs):
        handler = getattr(self, command)
        return handler(**kwargs)

    def build_job_submission(self, job=None, cfg=None, output_dir=None, 
                             **kwargs):
        """Build a job submission.

        Args:
            job (dict): job dict.
            cfg (dict): cfg dict.
            output_dir (str): dir in which to put submission files.

        Returns:
            submission_meta (dict): dict of submission_meta.

        """
        output_dir = output_dir or tempfile.mkdtemp()
        os.makedirs(output_dir, exist_ok=True)
        job_module = self.get_job_module(job=job, cfg=cfg)
        submission_meta_from_module = job_module.build_job_submission(
            job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
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

    def run_job_submission(self, submission_dir=None, **kwargs):
        """Run a job submission.

        Args:
            submission_dir (str): path to submission dir. This dir should
                contain the submission_meta file '{SUBMISSION_META_NAME}'
        """.format(SUBMISSION_META_NAME=self.SUBMISSION_META_NAME)
        submission_meta = self.read_submission_meta(dir_=submission_dir)
        job_module = self.get_job_module(job=submission_meta['job'],
                                         cfg=submission_meta['cfg'])
        return job_module.run_job_submission(submission_meta=submission_meta,
                                             **kwargs)

    def read_submission_meta(self, dir_=None):
        return _job_engine_utils.read_submission_meta(
            submission_dir=dir_, submission_meta_name=self.SUBMISSION_META_NAME)
