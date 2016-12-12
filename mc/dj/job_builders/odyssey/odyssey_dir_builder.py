import collections
import os
import tempfile
import textwrap


class OdysseyDirBuilder(object):
    def build_job(self, job_spec=None, output_dir=None, job_type_handlers=None):
        if not output_dir:
            output_dir = tempfile.mkdtemp(prefix='odyssey_dir.')
        dir_spec = self.generate_dir_spec(job_spec=job_spec,
                                          job_type_handlers=job_type_handlers)
        self.build_dir(dir_spec=dir_spec, output_dir=output_dir)
        return output_dir

    def generate_dir_spec(self, job_spec=None, job_type_handlers=None):
        dir_spec = self.generate_initial_dir_spec(job_spec=job_spec)
        job_type = job_spec['type']
        try:
            dir_spec_handler = job_type_handlers[job_type]
            dir_spec = dir_spec_handler(dir_spec=dir_spec, job_spec=job_spec)
        except KeyError as error:
            raise KeyError("No handler found for job_type '%s'." % job_type)
        return dir_spec

    def generate_initial_dir_spec(self, job_spec=None):
        dir_spec = collections.defaultdict(collections.OrderedDict)
        dir_spec['sbatch_params'] = self.generate_initial_sbatch_params(
            job_spec=job_spec)
        return dir_spec

    def generate_initial_sbatch_params(self, job_spec=None):
        initial_sbatch_params = collections.OrderedDict([
            ('nodes', 1),
            ('ntasks', 1),
            ('partition', 'my_partition')
        ])
        return initial_sbatch_params

    def build_dir(self, dir_spec=None, output_dir=None):
        self.write_job_script(dir_spec=dir_spec,
                              script_path=os.path.join(output_dir, 'job.sh'))

    def write_job_script(self, dir_spec=None, script_path=None):
        script_content = self.generate_job_script_content(dir_spec=dir_spec)
        with open(script_path, 'w') as f: f.write(script_content)

    def generate_job_script_content(self, dir_spec=None):
        content = textwrap.dedent(
            """
            #!/bin/bash
            {sbatch_section}

            {env_section}

            {module_section}

            {body_section}
            """
        ).format(
            sbatch_section=self.generate_sbatch_section(dir_spec=dir_spec),
            env_section=self.generate_env_section(dir_spec=dir_spec),
            module_section=self.generate_module_section(dir_spec=dir_spec),
            body_section=self.generate_body_section(dir_spec=dir_spec)
        ).lstrip()
        return content

    def generate_sbatch_section(self, dir_spec=None):
        sbatch_lines = ["#SBATCH --{key}={value}".format(key=key, value=value)
                        for key, value in dir_spec['sbatch_params'].items()]
        return "\n".join(sbatch_lines)


    def generate_env_section(self, dir_spec=None):
        env_lines = ["{key}={value}".format(key=key, value=value)
                     for key, value in dir_spec['env_vars'].items()]
        return "\n".join(env_lines)

    def generate_module_section(self, dir_spec=None):
        module_lines = ["module load {module}".format(module=module)
                        for module, include in dir_spec['modules'].items()
                        if include]
        return "\n".join(module_lines)

    def generate_body_section(self, dir_spec=None):
        return dir_spec.get('script_body', '')
