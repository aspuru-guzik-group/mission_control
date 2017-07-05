import json
import logging
import os
import tempfile

from . import constants
from .default_job_module_loader import DefaultJobModuleLoader


class JobModuleCommandDispatcher(object):
    SUBMISSION_META_NAME = constants.SUBMISSION_META_NAME

    def __init__(self, job_module_loader=None, logger=None):
        """
        Args:
            job_module_loader (JobModuleLoader): loader to use for loading
                job modules. Default: DefaultJobModuleLoader.
            logger: (logging.Logger): a logger.
        """
        self.logger = logger or logging
        self.job_module_loader = job_module_loader or DefaultJobModuleLoader()

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
        job_module = self.job_module_loader.load_job_module(job=job, cfg=cfg)
        submission_meta_from_module = job_module.build_job_submission(
            job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
        submission_meta = {'job': job, 'cfg': cfg, 'dir': output_dir,
                           **submission_meta_from_module}
        self._write_submission_meta(submission_meta=submission_meta,
                                    dir_=output_dir)
        return submission_meta

    def _write_submission_meta(self, submission_meta=None, dir_=None):
        submission_meta_path = os.path.join(dir_, self.SUBMISSION_META_NAME)
        with open(submission_meta_path, 'w') as f: json.dump(submission_meta, f)

    def run_job_submission(self, submission_dir=None,
                           submission_runner_cfg=None):
        """Run a job submission.

        Args:
            submission_dir (str): path to submission dir. This dir should
                contain the submission_meta file '{SUBMISSION_META_NAME}'
        """.format(SUBMISSION_META_NAME=self.SUBMISSION_META_NAME)
        submission_meta = self._read_submission_meta(dir_=submission_dir)
        job = submission_meta['job']
        cfg = self._get_run_submission_cfg(
            submission_runner_cfg=submission_runner_cfg, job=job)
        job_module = self.job_module_loader.load_job_module(job=job, cfg=cfg)
        return job_module.run_job_submission(submission_meta=submission_meta,
                                             cfg=cfg)

    def _read_submission_meta(self, dir_=None):
        """Read submission_meta from a submission dir."""
        submission_meta_path = os.path.join(dir_, self.SUBMISSION_META_NAME)
        with open(submission_meta_path) as f: return json.load(f)

    def _get_run_submission_cfg(self, submission_runner_cfg=None, job=None):
        get_cfg_fn = getattr(submission_runner_cfg, 'get_cfg_for_job', None)
        if get_cfg_fn: return get_cfg_fn(job=job)
        return {}
