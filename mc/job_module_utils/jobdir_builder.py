import os
import textwrap


class JobdirBuilder(object):
    """Class that builds jobdirs."""

    RUN_JOBDIR_CMD_VARNAME = 'MC_RUN_JOBDIR_CMD'
    ENTRYPOINT_NAME = 'entrypoint.sh'
    STD_LOG_FILE_NAMES = {
        log_name: (
            '{prefix}{log_name}'.format(
                prefix='MC.', log_name=log_name.upper())
        )
        for log_name in ['stdout', 'stderr', 'failure']
    }
    CHECKPOINT_FILE_NAMES = {
        checkpoint_name: (
            '{prefix}{checkpoint_name}'.format(
                prefix='MC.', checkpoint_name=checkpoint_name.upper())
        )
        for checkpoint_name in ['completed']
    }

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
                An alternative way to customize is to subclass and override
                :meth:`_decorate_job_spec`.
        """
        self.job = job or {}
        self.cfg = cfg or {}
        self.output_dir = output_dir
        self.entrypoint_path = os.path.join(self.output_dir,
                                            self.ENTRYPOINT_NAME)
        self.decorate_job_spec_fn = decorate_job_spec_fn

    @classmethod
    def build_jobdir(cls, job=None, cfg=None, output_dir=None, **kwargs):
        """Instantiates class and calls :meth:`_build_jobdir`."""
        instance = cls(job=job, cfg=cfg, output_dir=output_dir, **kwargs)
        return instance._build_jobdir(**kwargs)

    def _build_jobdir(self, **kwargs):
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
            {checkpoint_section}
            ${RUN_JOBDIR_CMD_VARNAME}
            """
        ).lstrip().format(
            checkpoint_section=self._generate_checkpoint_section(),
            RUN_JOBDIR_CMD_VARNAME=self.RUN_JOBDIR_CMD_VARNAME
        )

    def _generate_checkpoint_section(self):
        def generate_log_summary_cmd(log_name):
            return 'tail -n 50 %s' % self.STD_LOG_FILE_NAMES[log_name]

        return textwrap.dedent(
            '''
            START_DIR=$PWD
            output_checkpoint_files () {{
                PREV_RETURN_CODE=$?
                pushd $START_DIR > /dev/null
                if [ $PREV_RETURN_CODE -eq 0 ]; then
                    touch "{completed_checkpoint_name}"
                else
                    echo "{summary_stdout_cmd}:" >> "{failure_log_name}"
                    {summary_stdout_cmd} >> "{failure_log_name}"
                    echo "{summary_stderr_cmd}:" >> "{failure_log_name}"
                    {summary_stderr_cmd} >> "{failure_log_name}"
                    echo "{ls_cmd}:" >> "{failure_log_name}"
                    {ls_cmd} >> "{failure_log_name}"
                fi
                popd > /dev/null
            }}
            trap "output_checkpoint_files" EXIT
            '''
        ).strip().format(
            completed_checkpoint_name=self.CHECKPOINT_FILE_NAMES['completed'],
            failure_log_name=self.STD_LOG_FILE_NAMES['failure'],
            summary_stdout_cmd=generate_log_summary_cmd(log_name='stdout'),
            summary_stderr_cmd=generate_log_summary_cmd(log_name='stderr'),
            ls_cmd='ls -1'
        )

    def _generate_job_spec(self):
        job_spec = self._generate_base_job_spec()
        job_spec = self._decorate_job_spec(job_spec=job_spec)
        return job_spec

    def _generate_base_job_spec(self):
        return {
            'cfg_specs': {
                self.RUN_JOBDIR_CMD_VARNAME: {'required': True}
            },
            'dir': self.output_dir,
            'entrypoint': './' + self.ENTRYPOINT_NAME,
            'std_log_file_names': self.STD_LOG_FILE_NAMES
        }

    def _decorate_job_spec(self, job_spec=None):
        """Subclasses of this class can overrde this method to define custom
        logic.

        See the decorate_job_spec_fn parameter of `:meth:__init__` for other
        ways to customize that don't involve subclassing.

        The default behavior is to call the decorate_job_spec_fn provided to
        the constructor if it was provided. Otherwise it returns the job_spec
        unaltered.
        """
        if self.decorate_job_spec_fn:
            job_spec = self.decorate_job_spec_fn(job_spec=job_spec,
                                                 job=self.job,
                                                 cfg=self.cfg)
        return job_spec

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
