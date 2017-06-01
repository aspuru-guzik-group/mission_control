import json
import os
import textwrap

from .. import constants as confgen_constants


class WorkdirBuilder(object):
    def __init__(self, workdir=None, workdir_params=None, entrypoint_name=None):
        self.workdir = workdir
        self.workdir_params = workdir_params
        self.entrypoint_name = entrypoint_name or 'entrypoint.sh'
        self.infile_name = confgen_constants.CONFGEN_IN_FILENAME
        self.outdir_name = 'outputs'

    def build_workdir(self):
        self.ensure_dir(dir=self.workdir)
        for component in ['entrypoint', 'infile']:
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
            mkdir {outdir}
            $CONFGEN_EXE --infile={infile} --outdir={outdir}
            '''
        ).lstrip().format(infile=self.infile_name, outdir=self.outdir_name)
        return entrypoint_content

    def generate_infile_content(self):
        confgen_params = self.workdir_params['confgen_params']
        return json.dumps(confgen_params)
