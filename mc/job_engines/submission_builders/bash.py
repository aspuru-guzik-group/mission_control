import collections
import os
import textwrap
import tempfile

from jinja2 import Template


class BashSubmissionBuilder(object):

    checkpoint_files = {
        'completed': '_JOB__COMPLETED',
        'failed': '_JOB__FAILED',
    }

    entrypoint_name = 'job.sh'

    io_dirs = {'inputs': 'inputs', 'outputs': 'outputs'}

    std_log_filenames = {
        'stdout': '_JOB.stdout',
        'stderr': '_JOB.stderr',
        'failure': checkpoint_files['failed'],
    }

    def build_submission(self, submission_spec=None, cfg=None, output_dir=None):
        self.spec = submission_spec
        self.cfg = cfg
        self.output_dir = output_dir or tempfile.mkdtemp(prefix='bash.sub.')
        self.std_log_files = {
            logname: os.path.join(self.output_dir, filename)
            for logname, filename in self.std_log_filenames.items()
        }
        self.ensure_dir(output_dir)
        self.setup_io_dirs()
        self.write_entrypoint()
        self.write_templates()
        submission_meta = {
            'dir': self.output_dir,
            'checkpoint_files': self.checkpoint_files,
            'entrypoint': self.entrypoint_name,
            'std_log_files': self.std_log_files,
            'inputs_dir': 'inputs',
            'outputs_dir': 'outputs',
        }
        return submission_meta

    def ensure_dir(self, dir_path=None):
        if not os.path.exists(dir_path): os.makedirs(dir_path)

    def setup_io_dirs(self):
        for dir_name, rel_path in self.io_dirs.items():
            os.makedirs(os.path.join(self.output_dir, rel_path), exist_ok=True)

    def write_entrypoint(self):
        entrypoint_path = os.path.join(self.output_dir, self.entrypoint_name)
        with open(entrypoint_path, 'w') as f:
            f.write(self.generate_entrypoint_content())
        os.chmod(entrypoint_path, 0o755)

    def generate_entrypoint_content(self):
        content = textwrap.dedent(
            """
            #!/bin/bash
            {checkpoint_file_section}

            {env_section}

            {body_section}
            """
        ).format(
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
            completed_checkpoint_file=self.checkpoint_files['completed'],
            failure_log_file=std_log_files['failure'],
            summary_stdout_cmd=generate_summary_cmd(std_log_files['stdout']),
            summary_stderr_cmd=generate_summary_cmd(std_log_files['stderr']),
            ls_cmd='ls -1'
        )

    def get_std_log_files(self):
        return {**self.std_log_files, **self.spec.get('std_log_files', {})}

    def generate_env_section(self):
        env_vars = self.squash_kvp_list(self.spec.get('env_vars', []))
        env_lines = ["{key}={value}".format(key=key, value=value)
                     for key, value in env_vars]
        return "\n".join(env_lines)

    def squash_kvp_list(self, kvp_list=None):
        squashed = collections.OrderedDict()
        for k, v in kvp_list: squashed[k] = v
        return squashed.items()

    def generate_body_section(self):
        job_engine_cfg = self.cfg.get('job_engine', {})
        body_section = textwrap.dedent(
            """
            {job_engine_preamble}
            {job_engine_exe}  run_submission --submission_dir={output_dir}
            """
        ).strip().format(
            job_engine_preamble=job_engine_cfg.get('entrypoint_preamble', ''),
            job_engine_exe=job_engine_cfg['job_engine_exe'],
            output_dir=self.output_dir
        )
        return body_section

    def write_templates(self):
        templates = self.spec.get('templates', {})
        if not templates.get('specs', None): return
        ctx = templates.get('ctx', {})
        for template_spec in templates.get('specs', []):
            self.write_template_spec(template_spec=template_spec, ctx=ctx)

    def write_template_spec(self, template_spec=None, ctx=None):
        target_path = os.path.join(self.output_dir, template_spec['target'])
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
