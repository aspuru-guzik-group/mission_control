import collections
import os
import tempfile
import textwrap
import unittest
from unittest.mock import call, patch, MagicMock

from mc.mc_utils import test_utils as mc_test_utils
from .. import qchem_job_module
from .. import chemthing_utils
from .test_base_command import BaseCommandBaseTestCase


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.chemthing = self.generate_chemthing()
        self.job = {
            'data': {
                'input': {
                    'chemthing': self.chemthing
                }
            }
        }
        self.cfg = {}
        self.ctx_dir = 'some_ctx_dir'
        self.output_dir = 'some_output_dir'
        self.engine = qchem_job_module.QChemJobEngine()

    def generate_chemthing(self):
        mols = chemthing_utils.generate_mols_w_conformers(
            smiles_list=['Cc1ccccc1'])
        conformer = mols[0].GetConformers()[0]
        chemthing = chemthing_utils.generate_chemthing_for_conformer(
            conformer=conformer)
        return chemthing

class ExecuteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.engine.generate_qchem_dir = MagicMock()
        self.engine.execute_dir = MagicMock()

    def _execute_job(self, **kwargs):
        default_kwargs = {
            'job': self.job,
            'cfg': self.cfg,
            'output_dir': self.output_dir,
            'ctx_dir': self.ctx_dir,
        }
        self.engine.execute_job(**{**default_kwargs, **kwargs})

    def test_generates_qchem_dir(self):
        self._execute_job()
        self.assertEqual(self.engine.generate_qchem_dir.call_args,
                         call(job=self.job, output_dir=self.output_dir))

    def test_executes_qchem_dir(self):
        self._execute_job()
        expected_dir_meta = self.engine.generate_qchem_dir.return_value
        self.assertEqual(self.engine.execute_dir.call_args,
                         call(dir_meta=expected_dir_meta))

class ExecuteDirTestCase(BaseTestCase):
    def test_executes_dir_entrypoint(self):
        mock_dir_meta = collections.defaultdict(MagicMock())
        with patch.object(qchem_job_module, 'subprocess') as mock_subprocess:
            self.engine.execute_dir(dir_meta=mock_dir_meta)
            expected_cmd = 'cd {dir} && bash {entrypoint}'.format(
                **mock_dir_meta)
            self.assertEqual(
                mock_subprocess.run.call_args,
                call(expected_cmd, shell=True, check=True)
            )

class GenerateQChemDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.output_dir = tempfile.mkdtemp()

    def test_generates_expected_dir_files(self):
        dir_meta = self.engine.generate_qchem_dir(job=self.job,
                                                  output_dir=self.output_dir)
        actual_dir = dir_meta['dir']
        expected_dir = self.generate_expected_dir()
        mc_test_utils.assert_dirs_equal(
            test=self,
            left=actual_dir,
            right=expected_dir
        )

    def generate_expected_dir(self):
        expected_dir = tempfile.mkdtemp()
        def _write_dir_file(filename=None, content=None):
            with open(os.path.join(expected_dir, filename), 'w') as f:
                f.write(content)
        _write_dir_file(
            filename='qchem.inp',
            content=textwrap.dedent(
                '''
                $molecule
                READ coords.xyz
                $end

                $rem
                MEM_TOTAL 12000
                MEM_STATIC 2500
                JOBTYPE OPT
                GEOM_OPT_MAX_CYCLES 200
                GEOM_OPT_PRINT 2

                SYMMETRY FALSE

                !BP86
                EXCHANGE B
                CORRELATION P86

                BASIS 6-31G*
                PURECART 1111

                SCF_GUESS CORE
                XC_GRID 1
                MAX_SCF_CYCLES 200

                PRINT_ORBITALS FALSE
                SCF_PRINT 1
                SCF_FINAL_PRINT 1
                POP_MULLIKEN 0

                $end
                '''
            )
        )
        _write_dir_file(
            filename='coords.xyz',
            content=chemthing_utils.from_cml(cml=self.chemthing['cml'],
                                             output_format='xyz')
        )
        _write_dir_file(
            filename='entrypoint.sh',
            content=textwrap.dedent(
                '''
                #!/bin/bash
                INPUT_FILE="qchem.inp"
                OUTPUT_FILE="qchem.out"
                qchem $INPUT_FILE $OUTPUT_FILE
                '''
            )
        )
        return expected_dir

    def test_returns_expected_dir_meta(self):
        dir_meta = self.engine.generate_qchem_dir(job=self.job,
                                                  output_dir=self.output_dir)
        expected_dir_meta = {
            'dir': self.output_dir,
            'entrypoint': 'entrypoint.sh',
        }
        self.assertEqual(dir_meta, expected_dir_meta)


class CommandTestCase(BaseCommandBaseTestCase):
    def generate_command(self): return qchem_job_module.Command()

    def test_calls_execute_job(self):
        with patch.object(qchem_job_module, 'QChemJobEngine') as MockEngine:
            self._execute_command()
            self.assertEqual(MockEngine.call_args, call())
            self.assertEqual(
                MockEngine.return_value.execute_job.call_args,
                call(job=self.job, cfg=self.cfg, ctx_dir=self.ctx_dir,
                     output_dir=self.output_dir)
            )
