import collections
import json
import os
import tempfile
import textwrap

from jinja2 import Template


class OdysseyJobDirBuilder(object):

    checkpoint_files = {
        'completed': 'ODYSSEY_JOB__COMPLETED',
        'failed': 'ODYSSEY_JOB__FAILED',
    }

    std_log_files = {
        'stdout': 'ODYSSEY_JOB.stdout',
        'stderr': 'ODYSSEY_JOB.stderr',
        'failure': checkpoint_files['failed'],
    }

    @classmethod
    def build_dir(cls, dir_spec=None, submission_dir=None,
                  submission_meta_file_name='submission.json'):
        if not submission_dir: submission_dir = tempfile.mkdtemp(prefix='ody.')
        if not os.path.exists(submission_dir): os.makedirs(submission_dir)
        cls.write_templates(dir_spec=dir_spec, submission_dir=submission_dir)
        cls.write_entrypoint(dir_spec=dir_spec, submission_dir=submission_dir)
        for io_dir in  ['inputs', 'outputs']:
            os.makedirs(os.path.join(submission_dir, io_dir), exist_ok=True)
        submission_meta = {
            'dir': submission_dir,
            'checkpoint_files': cls.checkpoint_files,
            'entrypoint': 'job.sh',
            'std_log_files': cls.std_log_files,
            'inputs_dir': 'inputs',
            'outputs_dir': 'outputs',
        }
        if submission_meta_file_name:
            cls.write_submission_meta(
                submission_dir=submission_dir,
                submission_meta_file_name=submission_meta_file_name,
                submission_meta=submission_meta)
        return submission_meta

    @classmethod
    def write_entrypoint(cls, dir_spec=None, submission_dir=None):
        script_path = os.path.join(submission_dir, 'job.sh')
        script_content = cls.generate_entrypoint_content(dir_spec=dir_spec)
        with open(script_path, 'w') as f: f.write(script_content)

    @classmethod
    def generate_entrypoint_content(cls, dir_spec=None):
        content = textwrap.dedent(
            """
            #!/bin/bash
            {sbatch_section}

            {status_file_section}

            {env_section}

            {module_section}

            {body_section}
            """
        ).format(
            sbatch_section=cls.generate_sbatch_section(dir_spec=dir_spec),
            status_file_section=cls.generate_status_file_section(
                dir_spec=dir_spec),
            env_section=cls.generate_env_section(dir_spec=dir_spec),
            module_section=cls.generate_module_section(dir_spec=dir_spec),
            body_section=cls.generate_body_section(dir_spec=dir_spec)
        ).lstrip()
        return content

    @classmethod
    def generate_sbatch_section(cls, dir_spec=None):
        std_log_files = cls.get_std_log_files(dir_spec=dir_spec)
        std_log_file_sbatch_params = [
            ('error', std_log_files['stderr']),
            ('output', std_log_files['stdout']),
        ]
        sbatch_params = dir_spec.get('sbatch', []) + std_log_file_sbatch_params
        squashed_sbatch_params = cls.squash_kvp_list(sbatch_params)
        sbatch_lines = ["#SBATCH --{key}={value}".format(key=key, value=value)
                        for key, value in squashed_sbatch_params]
        return "\n".join(sbatch_lines)

    @classmethod
    def get_std_log_files(cls, dir_spec=None):
        return {
            **cls.std_log_files,
            **dir_spec.get('std_log_files', {})
        }

    @classmethod
    def squash_kvp_list(cls, kvp_list=None):
        squashed = collections.OrderedDict()
        for k, v in kvp_list: squashed[k] = v
        return squashed.items()

    @classmethod
    def generate_status_file_section(cls, dir_spec=None):
        std_log_files = cls.get_std_log_files(dir_spec=dir_spec)

        def generate_summary_cmd(file_name):
            return 'tail -n 50 %s' % file_name

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
            completed_checkpoint_file=cls.checkpoint_files['completed'],
            failure_log_file=cls.std_log_files['failure'],
            summary_stdout_cmd=generate_summary_cmd(std_log_files['stdout']),
            summary_stderr_cmd=generate_summary_cmd(std_log_files['stderr']),
            ls_cmd='ls -1'
        )

    @classmethod
    def generate_env_section(cls, dir_spec=None):
        env_vars = cls.squash_kvp_list(dir_spec.get('env_vars', []))
        env_lines = ["{key}={value}".format(key=key, value=value)
                     for key, value in env_vars]
        return "\n".join(env_lines)

    @classmethod
    def generate_module_section(cls, dir_spec=None):
        modules = cls.prepare_modules(modules=dir_spec.get('modules'))
        module_lines = ["module load {module}".format(module=module)
                        for module, include in modules if include]
        return "\n".join(module_lines)

    @classmethod
    def prepare_modules(cls, modules=None):
        prepared_modules = []
        for module in (modules or []):
            if isinstance(module, str): module = [module, True]
            prepared_modules.append(module)
        prepared_modules = cls.squash_kvp_list(prepared_modules)
        return prepared_modules

    @classmethod
    def generate_body_section(cls, dir_spec=None):
        return dir_spec.get('entrypoint_body', '')

    @classmethod
    def write_templates(cls, dir_spec=None, submission_dir=None):
        templates = dir_spec.get('templates', {})
        if not templates.get('specs', None): return
        ctx = templates.get('ctx', {})
        for template_spec in templates.get('specs', []):
            cls.write_template_spec(template_spec=template_spec,
                                    submission_dir=submission_dir,
                                    ctx=ctx)

    @classmethod
    def write_template_spec(cls, template_spec=None, submission_dir=None,
                            ctx=None):
        target_path = os.path.join(submission_dir, template_spec['target'])
        cls.ensure_dir(dir_path=os.path.dirname(target_path))
        rendered = cls.render_template_spec(template_spec=template_spec,
                                            ctx=ctx)
        with open(target_path, 'w') as f: f.write(rendered + "\n")

    @classmethod
    def render_template_spec(cls, template_spec=None, ctx=None):
        if template_spec.get('no_render', None):
            rendered = template_spec['content']
        else:
            template = Template(template_spec['content'])
            template.environment.keep_trailing_newline = True
            merged_ctx = {**(ctx or {}), **template_spec.get('ctx', {})}
            rendered = template.render(merged_ctx)
        return rendered

    @classmethod
    def ensure_dir(cls, dir_path=None):
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    @classmethod
    def write_submission_meta(cls, submission_dir=None,
                              submission_meta_file_name=None,
                              submission_meta=None):
        meta_file_path = os.path.join(submission_dir, submission_meta_file_name)
        pruned_submission_meta = {k: v for k, v in submission_meta.items()
                                  if k not in ['dir']}
        open(meta_file_path, 'w').write(json.dumps(pruned_submission_meta))
