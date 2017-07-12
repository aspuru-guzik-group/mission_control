import os
import textwrap


class JobdirBuilder(object):
    """Class that builds jobdirs."""

    RUN_JOBDIR_CMD_VARNAME = 'MC_RUN_JOBDIR_CMD'
    ENTRYPOINT_NAME = 'entrypoint.sh'

    def __init__(self, job=None, cfg=None, output_dir=None,
                 decorate_job_spec_fn=None, **kwargs):
        """
        Args:
            job (dict): job dict.
            cfg (dict): cfg dict.
            output_dir (str): directory in which to put job files.
            decorate_job_spec_fn (fn): a function that decorates a job_spec.
                This fn is how users of this class can customize job_specs.
                See :meth:`_default_decorate_job_spec_fn` for more info.
        """
        self.job = job or {}
        self.cfg = cfg or {}
        self.output_dir = output_dir
        self.entrypoint_path = os.path.join(self.output_dir,
                                            self.ENTRYPOINT_NAME)
        self.decorate_job_spec_fn = decorate_job_spec_fn

    def build_jobdir(self, *args, **kwargs):
        """Build a jobdir.

        Returns:
            submission_meta (dict): dict of submission metadata.
        """
        self._write_entrypoint()
        job_spec = self._generate_job_spec()
        return job_spec

    def _write_entrypoint(self):
        with open(self.entrypoint_path, 'w') as f:
            f.write(self._generate_entrypoint_content())
        os.chmod(self.entrypoint_path, 0o755)

    def _generate_entrypoint_content(self):
        return textwrap.dedent(
            """
            #!/bin/bash
            ${RUN_JOBDIR_CMD_VARNAME}
            """
        ).lstrip().format(RUN_JOBDIR_CMD_VARNAME=self.RUN_JOBDIR_CMD_VARNAME)

    def _generate_job_spec(self):
        job_spec = self._generate_base_job_spec()
        if self.decorate_job_spec_fn:
            job_spec = self.decorate_job_spec_fn(job_spec=job_spec,
                                                 job=self.job,
                                                 cfg=self.cfg)
        return job_spec

    def _generate_base_job_spec(self):
        return {
            'cfg_specs': {
                self.RUN_JOBDIR_CMD_VARNAME: {'required': True}
            },
            'dir': self.output_dir,
            'entrypoint': self.ENTRYPOINT_NAME,
        }

    def _default_decorate_job_spec_fn(self, job_spec=None, job=None, cfg=None):
        """Stub method to illustrate signature for decorate_job_spec_fn

        Args:
            job_spec (dict): the job_spec dict to decorate.
            job (dict): an mc job dict
            cfg (dict): a cfg dict

        Returns:
            decorated_job_spec (dict): the decorated job_spec.
        """
        return {**job_spec}
