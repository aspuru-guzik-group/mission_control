import collections
import os
import textwrap
import tempfile

from jinja2 import Template


class BashSubmissionBuilder(object):
    def __init__(self, *args, checkpoint_file_names=None,
                 entrypoint_file_name=None, io_dir_names=None,
                 std_log_file_names=None, **kwargs):
        self.checkpoint_file_names = {
            'completed': 'MC_JOB__COMPLETED',
            'failed': 'MC_JOB__FAILED',
            **(checkpoint_file_names or {})
        }
        self.entrypoint_file_name = entrypoint_file_name or 'job.sh'
        self.io_dir_names = {
            'inputs': 'inputs',
            'outputs': 'outputs',
            **(io_dir_names or {})
        }
        self.std_log_file_names = {
            'stdout': 'MC_JOB.stdout',
            'stderr': 'MC_JOB.stderr',
            'failure': self.checkpoint_file_names['failed'],
            **(std_log_file_names or {})
        }

    def build_submission(self, *args, submission_spec=None, cfg=None,
                         output_dir=None, **kwargs):
        self._spec = submission_spec or {}
        self._cfg = cfg
        self._output_dir = output_dir or tempfile.mkdtemp(prefix='bash.sub.')
        self.ensure_dir(self._output_dir)
        self.setup_io_dirs()
        self.write_entrypoint()
        self.write_templates()
        submission_meta = {
            'dir': self._output_dir,
            'checkpoint_file_names': self.checkpoint_file_names,
            'entrypoint_file_name': self.entrypoint_file_name,
            'std_log_file_names': self.std_log_file_names,
            'io_dir_names': self.io_dir_names,
        }
        return submission_meta

    def ensure_dir(self, dir_path=None):
        if not os.path.exists(dir_path): os.makedirs(dir_path)

    def setup_io_dirs(self):
        for dir_name, rel_path in self.io_dir_names.items():
            os.makedirs(os.path.join(self._output_dir, rel_path), exist_ok=True)

    def write_entrypoint(self):
        entrypoint_path = os.path.join(self._output_dir,
                                       self.entrypoint_file_name)
        with open(entrypoint_path, 'w') as f:
            f.write(self.generate_entrypoint_content())
        os.chmod(entrypoint_path, 0o755)

    def generate_entrypoint_content(self):
        content = textwrap.dedent(
            """
            {shebang_line}
            {header}
            {checkpoint_file_section}

            {env_section}

            {body_section}
            """
        ).format(
            shebang_line=(self._spec.get('shebang_line') or '#!/bin/bash'),
            header=(self._spec.get('header') or ''),
            checkpoint_file_section=self.generate_checkpoint_file_section(),
            env_section=self.generate_env_section(),
            body_section=self.generate_body_section()
        ).lstrip()
        return content

    def generate_checkpoint_file_section(self):
        std_log_files = self.get_std_log_files()
        def generate_summary_cmd(file_name): return 'tail -n 50 %s' % file_name
        return textwrap.dedent(
            '''
            START_DIR=$PWD
            output_status_file () {{
                PREV_RETURN_CODE=$?
                pushd $START_DIR > /dev/null
                if [ $PREV_RETURN_CODE -eq 0 ]; then
                    touch {completed_checkpoint_file}
                else
                    touch {failure_log_file}
                    echo "{summary_stdout_cmd}:" >> {failure_log_file}
                    {summary_stdout_cmd} >> {failure_log_file}
                    echo "{summary_stderr_cmd}:" >> {failure_log_file}
                    {summary_stderr_cmd} >> {failure_log_file}
                    echo "{ls_cmd}:" >> {failure_log_file}
                    {ls_cmd} >> {failure_log_file}
                fi
                popd > /dev/null
            }}
            trap "output_status_file" EXIT
            '''
        ).strip().format(
            completed_checkpoint_file=self.checkpoint_file_names['completed'],
            failure_log_file=std_log_files['failure'],
            summary_stdout_cmd=generate_summary_cmd(std_log_files['stdout']),
            summary_stderr_cmd=generate_summary_cmd(std_log_files['stderr']),
            ls_cmd='ls -1'
        )

    def get_std_log_files(self):
        return {
            **{
                log_name: os.path.join(self._output_dir, log_file_name)
                for log_name, log_file_name in self.std_log_file_names.items()
            },
            **self._spec.get('std_log_files', {})
        }

    def generate_env_section(self):
        env_vars = self.squash_kvp_list(self._spec.get('env_vars', []))
        env_lines = ["export {key}={value}".format(key=key, value=value)
                     for key, value in env_vars]
        return "\n".join(env_lines)

    def squash_kvp_list(self, kvp_list=None):
        squashed = collections.OrderedDict()
        for k, v in kvp_list: squashed[k] = v
        return squashed.items()

    def generate_body_section(self):
        job_engine_cfg = self._cfg.get('job_engine', {})
        body_section = textwrap.dedent(
            """
            {job_engine_preamble}
            {job_engine_exe} run_submission --submission_dir={output_dir}
            """
        ).strip().format(
            job_engine_preamble=job_engine_cfg.get('entrypoint_preamble', ''),
            job_engine_exe=job_engine_cfg['job_engine_exe'],
            output_dir=self._output_dir
        )
        return body_section

    def write_templates(self):
        templates = self._spec.get('templates', {})
        if not templates.get('specs', None): return
        ctx = templates.get('ctx', {})
        for template_spec in templates.get('specs', []):
            self.write_template_spec(template_spec=template_spec, ctx=ctx)

    def write_template_spec(self, template_spec=None, ctx=None):
        target_path = os.path.join(self._output_dir, template_spec['target'])
        self.ensure_dir(dir_path=os.path.dirname(target_path))
        rendered = self.render_template_spec(
            template_spec=template_spec, ctx=ctx)
        with open(target_path, 'w') as f: f.write(rendered + "\n")

    def render_template_spec(self, template_spec=None, ctx=None):
        if template_spec.get('no_render', None):
            rendered = template_spec['content']
        else:
            template = Template(template_spec['content'])
            template.environment.keep_trailing_newline = True
            merged_ctx = {**(ctx or {}), **template_spec.get('ctx', {})}
            rendered = template.render(merged_ctx)
        return rendered
