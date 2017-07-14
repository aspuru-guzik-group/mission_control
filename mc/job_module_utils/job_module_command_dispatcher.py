import json
import logging
import os
import tempfile

from . import constants
from . import utils
from .default_job_module_loader import DefaultJobModuleLoader


class JobModuleCommandDispatcher(object):
    def __init__(self, load_job_module_fn=None, logger=None):
        """
        Args:
            load_job_module_fn (fn): A fn to use for loading job modules. See
                :meth:`_detault_load_job_module_fn for more info.
            logger: (logging.Logger): a logger.
        """
        self.logger = logger or logging
        self.load_job_module_fn = load_job_module_fn or \
                self._default_load_job_module_fn

    def _default_load_job_module_fn(self, job=None, cfg=None):
        """Load a job module.

        Args:
            job (dict): job dict
            cfg (dict): cfg dict

        Returns:
            job_module (job_module): loaded module
        """
        return DefaultJobModuleLoader().load_job_module(job=job, cfg=cfg)

    def build_jobdir(self, job=None, cfg=None, output_dir=None, **kwargs):
        """Build a jobdir.

        Args:
            job (dict): job dict.
            cfg (dict): cfg dict.
            output_dir (str): dir in which to put job files.

        Returns:
            job_spec (jobman job spec): a dict that contains (1) the built
                jobdir path, and (2) parameters for running it.
        """
        output_dir = output_dir or tempfile.mkdtemp()
        os.makedirs(output_dir, exist_ok=True)
        job_module = self.load_job_module_fn(job=job, cfg=cfg)
        job_spec = job_module.build_jobdir(
            job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
        job_spec = {'job': job, 'cfg': cfg, 'dir': output_dir, **job_spec}
        self._write_job_meta_files(dir_=output_dir, job_spec=job_spec, job=job)
        return job_spec

    def _write_job_meta_files(self, dir_=None, **kwargs):
        for key, meta_file_name in constants.JOB_META_FILE_NAMES.items():
            path = os.path.join(dir_, meta_file_name)
            with open(path, 'w') as f: json.dump(kwargs[key], f)
        
    def run_jobdir(self, jobdir=None, cfg=None):
        """Run a jobdir.

        Args:
            jobdir (str): path to jobdir. This dir should
                contain: (1) the jobman job_spec file
                '{job_spec}', and (2) the mc job file '{job}'.
        """.format(**constants.JOB_META_FILE_NAMES)
        job_metas = {
            key: self._load_from_jobdir_file(
                jobdir=jobdir, file_name=meta_file_name)
            for key, meta_file_name in constants.JOB_META_FILE_NAMES.items()
        }
        job_module = self._load_job_module(job=job_metas['job'], cfg=None)
        return job_module.run_jobdir(**job_metas, cfg=cfg)

    def _load_job_module(self, job=None, cfg=None):
        return self.load_job_module_fn(job=job, cfg=cfg)

    def _load_from_jobdir_file(self, jobdir=None, file_name=None):
        return utils.load_from_jobdir_file(jobdir=jobdir, file_name=file_name)
