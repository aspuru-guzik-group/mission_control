import collections
import os
import tempfile
import textwrap

from jinja2 import Template


class OdysseyJobDirBuilder(object):
    @classmethod
    def build_dir(cls, dir_spec=None, output_dir=None):
        if not output_dir: output_dir = tempfile.mkdtemp(prefix='odyssey.')
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        cls.write_job_script(dir_spec=dir_spec, output_dir=output_dir)
        cls.write_templates(dir_spec=dir_spec, output_dir=output_dir)
        return output_dir

    @classmethod
    def write_job_script(cls, dir_spec=None, output_dir=None):
        script_path = os.path.join(output_dir, 'job.sh')
        script_content = cls.generate_job_script_content(dir_spec=dir_spec)
        with open(script_path, 'w') as f: f.write(script_content)

    @classmethod
    def generate_job_script_content(cls, dir_spec=None):
        content = textwrap.dedent(
            """
            #!/bin/bash
            {sbatch_section}

            {env_section}

            {module_section}

            {body_section}
            """
        ).format(
            sbatch_section=cls.generate_sbatch_section(dir_spec=dir_spec),
            env_section=cls.generate_env_section(dir_spec=dir_spec),
            module_section=cls.generate_module_section(dir_spec=dir_spec),
            body_section=cls.generate_body_section(dir_spec=dir_spec)
        ).lstrip()
        return content

    @classmethod
    def generate_sbatch_section(cls, dir_spec=None):
        sbatch_params = cls.flatten_kvp_list(dir_spec.get('sbatch', []))
        sbatch_lines = ["#SBATCH --{key}={value}".format(key=key, value=value)
                        for key, value in sbatch_params]
        return "\n".join(sbatch_lines)

    @classmethod
    def flatten_kvp_list(cls, kvp_list=None):
        flattened = collections.OrderedDict()
        for k, v in kvp_list: flattened[k] = v
        return flattened.items()

    @classmethod
    def generate_env_section(cls, dir_spec=None):
        env_vars = cls.flatten_kvp_list(dir_spec.get('env_vars', []))
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
        prepared_modules = cls.flatten_kvp_list(prepared_modules)
        return prepared_modules

    @classmethod
    def generate_body_section(cls, dir_spec=None):
        return dir_spec.get('job_script_body', '')

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
