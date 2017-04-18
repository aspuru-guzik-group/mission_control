import io
import json
import os
import tempfile
import textwrap
import unittest
from unittest.mock import call, MagicMock, patch

from ... import constants as qchem_constants
from .. import qchem_computation_parser

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.fixtures = Fixtures()
        self.parser = qchem_computation_parser.QChemComputationParser()

class Fixtures(object):
    def generate_abbreviated_qchem_header(self):
        return textwrap.dedent(
            '''
            Process 0 of 1 is on holy2a21202.rc.fas.harvard.edu - thread support 0
            initial socket setup ...start
            initial socket setup ...done 
            now start server 0 ... 
                              Welcome to Q-Chem
                 A Quantum Leap Into The Future Of Chemistry


             Q-Chem 4.2, Q-Chem, Inc., Pleasanton, CA (2014).

             Y. Shao,  Z. Gan,  E. Epifanovsky,  A. T. B. Gilbert,  M. Wormit,  
             J. M. Herbert,  A. I. Krylov,  P. M. W. Gill,  M. Head-Gordon

             Contributors to earlier versions of Q-Chem not listed above: 
             R. D. Adamson,  J. Baker,  E. F. C. Byrd,  A. K. Chakraborty,  C.-L. Cheng,  
             H. Dachsel,  R. J. Doerksen,  G. Hawkins,  A. Heyden,  S. Hirata,  

             Q-Chem 4.2.0 for Intel X86 EM64T Linux

             Parts of Q-Chem use Armadillo 4.100.2 (Dirt Cruiser).
             http://arma.sourceforge.net/

             Q-Chem begins on Tue Jul 19 15:15:57 2016  

            holy2a21202.rc.fas.harvard.edu
            Host: 0

                 Scratch files written to /scratch/qchem25325//
             Jun314 |n|home06|pkrastev|src|QCHEM|4.2|qchem_trunk|trunk|libmdc|qstr_info.C -1
             Parallel job on  1  processors
            Finally everything over in PARseQInput

            --------------------------------------------------------------
            User input:
            --------------------------------------------------------------
            $molecule
               READ coord.xyz
            $end

            $rem
            JOBTYPE OPT
            EXCHANGE B
            CORRELATION P86
            BASIS 6-31G*
            $end
            --------------------------------------------------------------
            ''')

    def generate_completed_qchem_job_dir(self):
        completed_qchem_dir = tempfile.mkdtemp()
        outputs_dir = os.path.join(completed_qchem_dir,
                                   self.completed_qchem_dir_outputs_name)
        os.makedirs(outputs_dir)
        self.generate_submission_file(parent_dir=completed_qchem_dir,
                                      outputs_dir=outputs_dir)
        self.generate_completed_qchem_work_dir(outputs_dir=outputs_dir)
        return completed_qchem_dir

    def generate_submission_file(self, parent_dir=None, outputs_dir=None):
        submission = {'outputs_dir': outputs_dir}
        with open(os.path.join(parent_dir, 'submission.json'), 'w') as f:
            json.dump(submission, f)

    def generate_completed_qchem_work_dir(self, outputs_dir=None):
        outputs_dir = outputs_dir or tempfile.mkdtemp()
        work_dir = os.path.join(outputs_dir, qchem_constants.QCHEM_OUTPUTS_KEY)
        os.makedirs(work_dir)
        qchem_output_file_path = os.path.join(
            work_dir, qchem_constants.QCHEM_OUTPUT_FILE_NAME)
        open(qchem_output_file_path, 'w').write(
            self.generate_abbreviated_qchem_header())
        return work_dir

class ExtractComputationMetaTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_dir = MagicMock()
        for method_name in ['get_qchem_output_file',
                            'extract_command_meta', 
                            'extract_execution_meta']:
            setattr(self.parser, method_name, MagicMock())
        self.computation_meta = self.parser.extract_computation_meta(
            job_dir=self.job_dir)

    def test_has_command_meta(self):
        expected_qchem_output_file = \
                self.parser.get_qchem_output_file.return_value
        self.assertEqual(
            self.parser.extract_command_meta.call_args,
            call(qchem_output_file=expected_qchem_output_file))
        self.assertEqual(self.computation_meta['command_meta'],
                         self.parser.extract_command_meta.return_value)

    def test_has_execution_meta(self):
        self.assertEqual(self.parser.extract_execution_meta.call_args,
                         call(job_dir=self.job_dir))
        self.assertEqual(self.computation_meta['execution_meta'],
                         self.parser.extract_execution_meta.return_value)

class GetQChemOutputFileTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_dir = MagicMock()
        self.parser.get_completed_qchem_work_dir = MagicMock()
        
    @patch.object(qchem_computation_parser, 'open')
    def test_returns_handle_for_qchem_output_file(self, mock_open):
        qchem_output_file = self.parser.get_qchem_output_file(
            job_dir=self.job_dir)
        expected_work_dir = \
                self.parser.get_completed_qchem_work_dir.return_value
        expected_output_file_path = os.path.join(
            expected_work_dir, qchem_constants.QCHEM_OUTPUT_FILE_NAME)
        self.assertEqual(qchem_output_file,
                         mock_open(expected_output_file_path))

class GetCompletedQChemWorkDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_dir = MagicMock()
        self.parser.parse_submission = MagicMock()
        
    def test_returns_path_to_work_dir(self):
        work_dir = self.parser.get_completed_qchem_work_dir(
            job_dir=self.job_dir)
        self.assertEqual(self.parser.parse_submission.call_args,
                         call(job_dir=self.job_dir))
        expected_submission = self.parser.parse_submission.return_value
        expected_work_dir = os.path.join(self.job_dir,
                                         expected_submission['outputs_dir'], 
                                         qchem_constants.QCHEM_OUTPUTS_KEY)
        self.assertEqual(work_dir, expected_work_dir)

class ParseSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_dir = MagicMock()
        
    @patch.object(qchem_computation_parser, 'open')
    @patch.object(qchem_computation_parser, 'json')
    def test_returns_path_to_work_dir(self, mock_json, mock_open):
        submission = self.parser.parse_submission(job_dir=self.job_dir)
        self.assertEqual(mock_json.load.call_args, call(mock_open.return_value))
        expected_submission_path = os.path.join(self.job_dir, 'submission.json')
        self.assertEqual(mock_open.call_args, call(expected_submission_path))
        self.assertEqual(submission, mock_json.load.return_value)

class ExtractCommandMetaTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.qchem_output_file = io.StringIO(
            self.fixtures.generate_abbreviated_qchem_header())
        self.command_meta = self.parser.extract_command_meta(
            qchem_output_file=self.qchem_output_file)

    def test_has_executable_meta(self):
        expected_executable_meta = {'name': 'Q-Chem', 'version': '4.2.0'}
        self.assertEqual(self.command_meta['executable'], 
                         expected_executable_meta)

    def test_has_qchem_params(self):
        expected_qchem_params = {
            'JOBTYPE': 'OPT',
            'EXCHANGE': 'B',
            'CORRELATION': 'P86',
            'BASIS': '6-31G*',
        }
        self.assertEqual(self.command_meta['qchem_params'], 
                         expected_qchem_params)

class ExtractExecutionMetaTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_dir = MagicMock()
        self.parser.parse_submission = MagicMock()

    @patch.object(qchem_computation_parser, 'open')
    @patch.object(qchem_computation_parser, 'json')
    def test_returns_execution_meta(self, mock_json, mock_open):
        execution_meta = self.parser.extract_execution_meta(
            job_dir=self.job_dir)
        self.assertEqual(mock_json.load.call_args, call(mock_open.return_value))
        expected_submission = self.parser.parse_submission.return_value
        expected_execution_meta_path = os.path.join(
            self.job_dir,
            expected_submission['outputs_dir'],
            expected_submission['execution_meta_name']
        )
        self.assertEqual(mock_open.call_args,
                         call(expected_execution_meta_path))
        self.assertEqual(execution_meta, mock_json.load.return_value)

class GenerateChemThingTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.computation_meta = MagicMock()
        self.parse_params = MagicMock()

    def _generate_chemthing(self):
        self.chemthing = self.parser.generate_chemthing(
            computation_meta=self.computation_meta,
            parse_params=self.parse_params)

    @patch.object(qchem_computation_parser, 'uuid4')
    def test_has_uuid(self, mock_uuid4):
        self._generate_chemthing()
        self.assertEqual(self.chemthing['uuid'], str(mock_uuid4.return_value))

    def test_has_artifacts(self):
        self._generate_chemthing()
        self.assertEqual(self.chemthing['props']['a2g2:prop:artifacts'],
                         self.parse_params['artifacts'])

    def test_has_command_meta(self):
        self._generate_chemthing()
        self.assertEqual(self.chemthing['props']['a2g2:prop:command_meta'],
                         self.computation_meta['command_meta'])

    def test_has_execution_meta(self):
        self._generate_chemthing()
        self.assertEqual(self.chemthing['props']['a2g2:prop:execution_meta'],
                         self.computation_meta['execution_meta'])

    def test_has_precursors(self):
        self._generate_chemthing()
        self.assertEqual(self.chemthing['precursors'],
                         self.parse_params.get('precursors'))

    def test_has_ancestors(self):
        self._generate_chemthing()
        self.assertEqual(self.chemthing['ancestors'],
                         self.parse_params.get('ancestors'))

class WriteChemThing(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.chemthing = MagicMock()
        self.chemthings_bulk_path = MagicMock()

    @patch.object(qchem_computation_parser, 'open')
    @patch.object(qchem_computation_parser, 'json')
    def test_writes_chemthings_bulk_file(self, mock_json, mock_open):
        self.parser.write_chemthing_bulk_file(chemthings=self.chemthings,
                                              path=self.chemthings_bulk_path)
        self.assertEqual(mock_json.dump.call_args,
                         call(self.chemthing, mock_open.return_value))
        self.assertEqual(mock_open.call_args, call(self.path))
