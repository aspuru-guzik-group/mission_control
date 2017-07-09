import json
import logging
import os
import tempfile

from . import constants
from .default_job_module_loader import DefaultJobModuleLoader


class JobModuleCommandDispatcher(object):
    JOB_SPEC_NAME = constants.JOB_SPEC_NAME

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
            job_spec (dict): a spec that describes where the built dir is,
                and how to run it.
        """
        output_dir = output_dir or tempfile.mkdtemp()
        os.makedirs(output_dir, exist_ok=True)
        job_module = self.job_module_loader.load_job_module(job=job, cfg=cfg)
        job_spec = job_module.build_jobdir(
            job=job, cfg=cfg, output_dir=output_dir, **kwargs) or {}
        job_spec = {'job': job, 'cfg': cfg, 'dir': output_dir, **job_spec}
        self._write_job_spec(job_spec=job_spec, dir_=output_dir)
        return job_spec

    def _write_job_spec(self, job_spec=None, dir_=None):
        job_spec_path = os.path.join(dir_, self.JOB_SPEC_NAME)
        with open(job_spec_path, 'w') as f: json.dump(job_spec, f)

    def run_jobdir(self, jobdir=None):
        """Run a jobdir.

        Args:
            jobdir (str): path to jobdir. This dir should
                contain the job_spec file '{JOB_SPEC_NAME}'
             (str): path to jobdir. This dir should
                contain the job_spec file '{JOB_SPEC_NAME}'
        """.format(JOB_SPEC_NAME=self.JOB_SPEC_NAME)
        job_spec = self._read_job_spec(dir_=jobdir)
        job = job_spec['job']
        job_module = self.job_module_loader.load_job_module(job=job)
        return job_module.run_jobdir(job_spec=job_spec)

    def _read_job_spec(self, dir_=None):
        """Read job_spec from a submission dir."""
        job_spec_path = os.path.join(dir_, self.JOB_SPEC_NAME)
        with open(job_spec_path) as f: return json.load(f)
