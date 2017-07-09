import json
import logging
import os
import tempfile

from . import constants
from .default_job_module_loader import DefaultJobModuleLoader


class JobModuleCommandDispatcher(object):
    JOBDIR_META_NAME = constants.JOBDIR_META_NAME

    def __init__(self, job_module_loader=None, logger=None):
        """
        Args:
            job_module_loader (JobModuleLoader): loader to use for loading
                job modules. Default: DefaultJobModuleLoader.
            logger: (logging.Logger): a logger.
        """
        self.logger = logger or logging
        self.job_module_loader = job_module_loader or DefaultJobModuleLoader()

    def build_jobdir(self, job=None, cfg=None, output_dir=None, **kwargs):
        """Build a jobdir.

        Args:
            job (dict): job dict.
            cfg (dict): cfg dict.
            output_dir (str): dir in which to put job files.

        Returns:
            jobdir_meta (dict): dict of jobdir_meta.

        """
        output_dir = output_dir or tempfile.mkdtemp()
        os.makedirs(output_dir, exist_ok=True)
        job_module = self.job_module_loader.load_job_module(job=job, cfg=cfg)
        jobdir_meta_from_module = job_module.build_jobdir(
            job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
        jobdir_meta = {'job': job, 'cfg': cfg, 'dir': output_dir,
                           **jobdir_meta_from_module}
        self._write_jobdir_meta(jobdir_meta=jobdir_meta,
                                    dir_=output_dir)
        return jobdir_meta

    def _write_jobdir_meta(self, jobdir_meta=None, dir_=None):
        jobdir_meta_path = os.path.join(dir_, self.JOBDIR_META_NAME)
        with open(jobdir_meta_path, 'w') as f: json.dump(jobdir_meta, f)

    def run_jobdir(self, jobdir=None):
        """Run a jobdir.

        Args:
            jobdir (str): path to jobdir. This dir should
                contain the jobdir_meta file '{JOBDIR_META_NAME}'
             (str): path to jobdir. This dir should
                contain the jobdir_meta file '{JOBDIR_META_NAME}'
        """.format(JOBDIR_META_NAME=self.JOBDIR_META_NAME)
        jobdir_meta = self._read_jobdir_meta(dir_=jobdir)
        job = jobdir_meta['job']
        job_module = self.job_module_loader.load_job_module(job=job)
        return job_module.run_jobdir(jobdir_meta=jobdir_meta)

    def _read_jobdir_meta(self, dir_=None):
        """Read jobdir_meta from a submission dir."""
        jobdir_meta_path = os.path.join(dir_, self.JOBDIR_META_NAME)
        with open(jobdir_meta_path) as f: return json.load(f)
