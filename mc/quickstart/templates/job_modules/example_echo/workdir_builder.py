import os
import textwrap

from . import cfg


class WorkdirBuilder(object):
    DEFAULT_ECHO_EXE = 'echo'

    class InvalidWorkdirParamsError(Exception): pass

    def __init__(self, workdir=None, workdir_params=None, entrypoint_name=None):
        self.workdir = workdir
        self.validate_workdir_params(workdir_params=workdir_params)
        self.workdir_params = workdir_params
        self.entrypoint_name = entrypoint_name or 'entrypoint.sh'
        self.input_file_name = cfg.INPUT_FILE_NAME
        self.output_file_name = cfg.OUTPUT_FILE_NAME

    def validate_workdir_params(self, workdir_params=None):
        try: assert workdir_params['message'] is not None
        except: raise self.InvalidWorkdirParamsError("Invalid message")

    def build_workdir(self):
        self.ensure_dir(dir_=self.workdir)
        for component in ['entrypoint', 'input_file']:
            path = os.path.join(self.workdir,
                                getattr(self, '%s_name' % component))
            content_generator = getattr(self, 'generate_%s_content' % component)
            with open(path, 'w') as f: f.write(content_generator())
            if component == 'entrypoint': os.chmod(path, 0o755)
        workdir_meta = {
            'dir': self.workdir,
            'entrypoint': self.entrypoint_name,
            'output_file_name': self.output_file_name,
        }
        return workdir_meta

    def ensure_dir(self, dir_=None):
        os.makedirs(dir_, exist_ok=True)

    def generate_entrypoint_content(self):
        entrypoint_content = textwrap.dedent(
            '''
            #!/bin/bash

            DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
            cd $DIR

            ${echo_cmd_output_key} "$(cat "{input_file_name}")" \\
                > "{output_file_name}"
            '''
        ).lstrip().format(
            echo_cmd_output_key=(
                cfg.cfg_specs['EXAMPLE_ECHO_ECHO_CMD']['output_key']),
            input_file_name=self.input_file_name,
            output_file_name=self.output_file_name
        )
        return entrypoint_content

    def generate_input_file_content(self): return self.workdir_params['message']
