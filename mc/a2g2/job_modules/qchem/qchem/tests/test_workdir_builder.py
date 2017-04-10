import os
import textwrap
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from mc.mc_utils import test_utils as mc_test_utils

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
        self.workdir_components = ['entrypoint', 'input_file']
        for component in self.workdir_components:
            content_generator_name = 'generate_%s_content' % component
            setattr(self.workdir_builder, content_generator_name,
                    MagicMock(return_value='%s_content' % component))

    def test_creates_expected_dir(self):
        mc_test_utils.assert_dirs_equal(
            test=self,
            left=self.workdir_builder.build_workdir()['dir'],
            right=self.generate_expected_dir()
        )

    def generate_expected_dir(self):
        expected_dir = tempfile.mkdtemp()
        for component in self.workdir_components:
            path = os.path.join(
                expected_dir,
                getattr(self.workdir_builder, '%s_name' % component))
            content_generator = getattr(
                self.workdir_builder, 'generate_%s_content' % component)
            open(path, 'w').write(content_generator())
        return expected_dir

    def test_returns_workdir_meta(self):
        workdir_meta = self.workdir_builder.build_workdir()
        expected_workdir_meta = {
            'dir': self.workdir_builder.workdir,
            'entrypoint': self.workdir_builder.entrypoint_name,
        }
        self.assertEqual(workdir_meta, expected_workdir_meta)

class GenerateEntrypointContentTestCase(BaseTestCase):
    def test_generates_expected_entrypoint_content(self):
        expected_entrypoint_content = textwrap.dedent(
            '''
            #!/bin/bash
            DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
            cd $DIR
            $QCHEM_EXE {infile} {outfile}
            '''
        ).lstrip().format(infile=self.workdir_builder.input_file_name,
                          outfile=self.workdir_builder.output_file_name)
        self.assertEqual(self.workdir_builder.generate_entrypoint_content(),
                         expected_entrypoint_content)

class GenerateInputFileContentTestCase(BaseTestCase):
    @patch.object(workdir_builder, 'qchem_input_generator')
    def test_generates_expected_content(self, mock_input_generator):
        self.assertEqual(self.workdir_builder.generate_input_file_content(),
                         mock_input_generator.generate_qchem_input.return_value)
