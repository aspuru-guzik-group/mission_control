import collections
import os
import tempfile
import textwrap


class OdysseyJobDirBuilder(object):
    @classmethod
    def build_dir(cls, dir_spec=None, output_dir=None):
        if not output_dir: output_dir = tempfile.mkdtemp(prefix='odyssey.')
        cls.write_job_script(dir_spec=dir_spec,
                             script_path=os.path.join(output_dir, 'job.sh'))

    @classmethod
    def write_job_script(cls, dir_spec=None, script_path=None):
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
        modules = cls.flatten_kvp_list(dir_spec.get('modules', []))
        module_lines = ["module load {module}".format(module=module)
                        for module, include in modules if include]
        return "\n".join(module_lines)

    @classmethod
    def generate_body_section(cls, dir_spec=None):
        return dir_spec.get('job_script_body', '')
