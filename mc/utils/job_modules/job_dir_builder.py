import json
from pathlib import Path
import textwrap

from a2g2_v2.utils.dispatcher import Dispatcher

from . import constants
from . import utils


class JobDirBuilder(object):
    FILE_COMPONENTS = ['job_key', 'job_meta', 'job_spec', 'entrypoint']

    class WorkdirError(Exception):
        def __init__(self, msg=None, job_dict=None, **kwargs):
            if not msg:
                msg = ''
            else:
                msg += ' '
            msg += (
                "Could not build work_dir; job_dict was '{job_dict}'"
            ).format(job_dict=job_dict)
            super().__init__(msg, **kwargs)

    def __init__(self, dispatcher=None):
        self.dispatcher = dispatcher or Dispatcher()

    def build_job_dir(self, job_dict=None, output_dir=None):
        self.job_dict = job_dict
        self.output_dir = output_dir
        self._ensure_dir(self.output_dir)
        self.work_dir_meta = self._build_work_dir()
        for component in self.FILE_COMPONENTS:
            self._build_file_component(component)
        return output_dir

    def _ensure_dir(self, dir=None):
        Path(dir).mkdir(parents=True, exist_ok=True)

    def _build_work_dir(self):
        try:
            work_dir_meta = self.dispatcher.dispatch(
                module_name=self.job_dict['module'], command='build_work_dir',
                params=self.job_dict['params'],
                output_dir=utils.get_job_dir_component_path(
                    job_dir=self.output_dir, component_name='work_dir')
            )
            return work_dir_meta
        except Exception as exc:
            raise self.WorkdirError(job_dict=self.job_dict) from exc

    def _build_file_component(self, component=None):
        content_fn = getattr(self, '_generate_%s_content' % component)
        content = content_fn()
        component_path = utils.get_job_dir_component_path(
            job_dir=self.output_dir, component_name=component)
        with open(component_path, 'w') as f:
            f.write(content)

    def _generate_job_key_content(self): return self.job_dict.get('key')

    def _generate_job_meta_content(self): return json.dumps(self.job_dict)

    def _generate_job_spec_content(self):
        return json.dumps(self._generate_job_spec())

    def _generate_job_spec(self):
        return {
            'entrypoint_name': constants.JOB_DIR_COMPONENT_PATHS['entrypoint'],
            'cfg_specs': self.work_dir_meta.get('cfg_specs', {})
        }

    def _generate_entrypoint_content(self):
        return textwrap.dedent(
            """
            #!/bin/bash
            DIR="{bash_get_src_dir}"
            pushd $DIR > /dev/null
            {checkpoint_section}
            {work_dir}/{work_dir_entrypoint_name} \\
                > >(tee -a {stdout_log}) \\
                2> >(tee -a {stderr_log} >&2)
            popd > /dev/null
            """
        ).lstrip().format(
            bash_get_src_dir='$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)',
            checkpoint_section=self._generate_checkpoint_section(),
            work_dir=constants.JOB_DIR_COMPONENT_PATHS['work_dir'],
            work_dir_entrypoint_name=self.work_dir_meta['entrypoint_name'],
            stdout_log=constants.JOB_DIR_COMPONENT_PATHS['stdout_log'],
            stderr_log=constants.JOB_DIR_COMPONENT_PATHS['stderr_log']
        )

    def _generate_checkpoint_section(self):
        def generate_log_summary_cmd(log_name):
            cmd = 'tail -n 50 %s' % constants.JOB_DIR_COMPONENT_PATHS[log_name]
            return cmd

        return textwrap.dedent(
            '''
            START_DIR=$PWD
            output_checkpoint_files () {{
                PREV_RETURN_CODE=$?
                pushd $START_DIR > /dev/null
                touch {executed_checkpoint}
                if [ $PREV_RETURN_CODE -eq 0 ]; then
                    touch "{completed_checkpoint}"
                else
                    echo "{summary_stdout_cmd}:" >> "{failed_checkpoint}"
                    {summary_stdout_cmd} >> "{failed_checkpoint}"
                    echo "{summary_stderr_cmd}:" >> "{failed_checkpoint}"
                    {summary_stderr_cmd} >> "{failed_checkpoint}"
                    echo "{ls_cmd}:" >> "{failed_checkpoint}"
                    {ls_cmd} >> "{failed_checkpoint}"
                fi
                popd > /dev/null
            }}
            trap "output_checkpoint_files" EXIT
            '''
        ).lstrip().format(
            executed_checkpoint=(
                constants.JOB_DIR_COMPONENT_PATHS['executed_checkpoint']),
            completed_checkpoint=(
                constants.JOB_DIR_COMPONENT_PATHS['completed_checkpoint']),
            failed_checkpoint=(
                constants.JOB_DIR_COMPONENT_PATHS['failed_checkpoint']),
            summary_stdout_cmd=generate_log_summary_cmd(log_name='stdout_log'),
            summary_stderr_cmd=generate_log_summary_cmd(log_name='stderr_log'),
            ls_cmd='ls -1',
        )
