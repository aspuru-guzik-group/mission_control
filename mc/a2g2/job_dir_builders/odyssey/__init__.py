import collections
import os
import tempfile
import textwrap

from jinja2 import Template


class OdysseyJobDirBuilder(object):

    checkpoint_files = {
        'completed': 'ODYSSEY_JOB__COMPLETED',
        'failed': 'ODYSSEY_JOB__FAILED',
    }

    output_files = {
        'stdout': 'ODYSSEY_JOB.stdout',
        'stderr': 'ODYSSEY_JOB.stderr',
    }

    @classmethod
    def build_dir(cls, dir_spec=None, output_dir=None):
        if not output_dir: output_dir = tempfile.mkdtemp(prefix='odyssey.')
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        cls.write_templates(dir_spec=dir_spec, output_dir=output_dir)
        cls.write_entrypoint(dir_spec=dir_spec, output_dir=output_dir)
        dir_meta = {
            'dir': output_dir,
            'checkpoint_files': cls.checkpoint_files,
            'entrypoint': 'job.sh',
            'output_files': cls.output_files,
        }
        return dir_meta

    @classmethod
    def write_entrypoint(cls, dir_spec=None, output_dir=None):
        script_path = os.path.join(output_dir, 'job.sh')
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
        output_files = cls.get_output_files(dir_spec=dir_spec)
        output_file_sbatch_params = [
            ('error', output_files['stderr']),
            ('output', output_files['stdout']),
        ]
        sbatch_params = dir_spec.get('sbatch', []) + output_file_sbatch_params
        squashed_sbatch_params = cls.squash_kvp_list(sbatch_params)
        sbatch_lines = ["#SBATCH --{key}={value}".format(key=key, value=value)
                        for key, value in squashed_sbatch_params]
        return "\n".join(sbatch_lines)

    @classmethod
    def get_output_files(cls, dir_spec=None):
        return {
            **cls.output_files,
            **dir_spec.get('output_files', {})
        }

    @classmethod
    def squash_kvp_list(cls, kvp_list=None):
        squashed = collections.OrderedDict()
        for k, v in kvp_list: squashed[k] = v
        return squashed.items()

    @classmethod
    def generate_status_file_section(cls, dir_spec=None):
        output_files = cls.get_output_files(dir_spec=dir_spec)

        def generate_tail_cmd(file_name):
            return 'tail -n 50 %s' % file_name

        return textwrap.dedent(
            '''
            START_DIR=$PWD
            output_status_file () {{
                PREV_RETURN_CODE=$?
                pushd $START_DIR
                if [ $PREV_RETURN_CODE -eq 0 ]; then
                    touch {completed_checkpoint_file}
                else
                    touch {failed_checkpoint_file}
                    echo "{tail_stdout_cmd}:" >> {failed_checkpoint_file}
                    {tail_stdout_cmd} >> {failed_checkpoint_file}
                    echo "{tail_stderr_cmd}:" >> {failed_checkpoint_file}
                    {tail_stderr_cmd} >> {failed_checkpoint_file}
                    echo "{ls_cmd}:" >> {failed_checkpoint_file}
                    {ls_cmd} >> {failed_checkpoint_file}
                fi
                popd
            }}
            trap "output_status_file" EXIT
            '''
        ).strip().format(
            completed_checkpoint_file=cls.checkpoint_files['completed'],
            failed_checkpoint_file=cls.checkpoint_files['failed'],
            tail_stdout_cmd=generate_tail_cmd(output_files['stdout']),
            tail_stderr_cmd=generate_tail_cmd(output_files['stderr']),
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
    def write_templates(cls, dir_spec=None, output_dir=None):
        templates = dir_spec.get('templates', {})
        if not templates.get('specs', None): return
        ctx = templates.get('ctx', {})
        for template_spec in templates.get('specs', []):
            cls.write_template_spec(template_spec=template_spec,
                                    output_dir=output_dir,
                                    ctx=ctx)

    @classmethod
    def write_template_spec(cls, template_spec=None, output_dir=None, ctx=None):
        target_path = os.path.join(output_dir, template_spec['target'])
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
