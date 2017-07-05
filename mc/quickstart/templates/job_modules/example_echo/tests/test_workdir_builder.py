import json
import os
import textwrap
import tempfile
import unittest
from unittest.mock import MagicMock

from a2g2_mc.utils import test_utils

from .. import workdir_builder


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.workdir = MagicMock()
        self.workdir_params = MagicMock()
        self.workdir_builder = workdir_builder.WorkdirBuilder(
            workdir=self.workdir, workdir_params=self.workdir_params)

class BuildWorkdirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.workdir_builder.workdir = tempfile.mkdtemp()
        self.workdir_components = ['entrypoint', 'infile']
        for component in self.workdir_components:
            content_generator_name = 'generate_%s_content' % component
            setattr(self.workdir_builder, content_generator_name,
                    MagicMock(return_value='%s_content' % component))

    def test_creates_expected_dir(self):
        test_utils.assert_dirs_equal(
            test=self,
            left=self.workdir_builder.build_workdir()['dir'],
            right=self.generate_expected_dir()
        )

    def test_returns_workdir_meta(self):
        workdir_meta = self.workdir_builder.build_workdir()
        expected_workdir_meta = {
            'dir': self.workdir_builder.workdir,
            'entrypoint': self.workdir_builder.entrypoint_name,
            'outdir_name': self.workdir_builder.outdir_name,
        }
        self.assertEqual(workdir_meta, expected_workdir_meta)

    def generate_expected_dir(self):
        expected_dir = tempfile.mkdtemp()
        for component in self.workdir_components:
            path = os.path.join(
                expected_dir,
                getattr(self.workdir_builder, '%s_name' % component))
            content_generator = getattr(
                self.workdir_builder, 'generate_%s_content' % component)
            with open(path, 'w') as f: f.write(content_generator())
        return expected_dir

class GenerateEntrypointContentTestCase(BaseTestCase):
    def test_generates_expected_entrypoint_content(self):
        expected_entrypoint_content = textwrap.dedent(
            '''
            #!/bin/bash
            DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
            cd $DIR
            mkdir {outdir}
            ${{CONFGEN_EXE:-{default_exe}}} \\
                --params_file={infile} --output_dir={outdir}
            '''
        ).lstrip().format(
            default_exe=self.workdir_builder.DEFAULT_EXE,
            infile=self.workdir_builder.infile_name,
            outdir=self.workdir_builder.outdir_name
        )
        self.assertEqual(self.workdir_builder.generate_entrypoint_content(),
                         expected_entrypoint_content)

class GenerateInfileContentTestCase(BaseTestCase):
    def test_generates_expected_infile_content(self):
        confgen_params = {'some': 'confgen_params'}
        self.workdir_builder.workdir_params = {
            'confgen_params': confgen_params
        }
        expected_infile_content = json.dumps(confgen_params)
        self.assertEqual(self.workdir_builder.generate_infile_content(),
                         expected_infile_content)
