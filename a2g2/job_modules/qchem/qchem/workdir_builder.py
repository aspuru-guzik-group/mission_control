import os
import textwrap

from .. import constants as qchem_constants
from . import qchem_input_generator


class WorkdirBuilder(object):
    def __init__(self, workdir=None, workdir_params=None, entrypoint_name=None):
        self.workdir = workdir
        self.workdir_params = workdir_params
        self.entrypoint_name = entrypoint_name or 'entrypoint.sh'
        self.input_file_name = qchem_constants.QCHEM_INPUT_FILE_NAME
        self.output_file_name = qchem_constants.QCHEM_OUTPUT_FILE_NAME

    def build_workdir(self):
        self.ensure_dir(dir=self.workdir)
        for component in ['entrypoint', 'input_file']:
            path = os.path.join(self.workdir,
                                getattr(self, '%s_name' % component))
            content_generator = getattr(self, 'generate_%s_content' % component)
            open(path, 'w').write(content_generator())
        workdir_meta = {
            'dir': self.workdir,
            'entrypoint': self.entrypoint_name,
        }
        return workdir_meta

    def ensure_dir(self, dir=None):
        os.makedirs(dir, exist_ok=True)

    def generate_entrypoint_content(self):
        entrypoint_content = textwrap.dedent(
            '''
            #!/bin/bash
            DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
            cd $DIR
            $QCHEM_EXE {infile} {outfile}
            '''
        ).lstrip().format(infile=self.input_file_name,
                          outfile=self.output_file_name)
        return entrypoint_content

    def generate_input_file_content(self):
        return qchem_input_generator.generate_qchem_input(
            params=self.workdir_params['qchem_params'])
