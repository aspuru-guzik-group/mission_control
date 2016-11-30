import glob
import os
import tempfile
from unittest.mock import call, patch, DEFAULT

from django.test import TestCase

from ..models import WorkflowRunner
from .. import utils

class SyncWorkflowRunnersTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, add_workflow_runner=DEFAULT,
                                    get_default_runner_dirs=DEFAULT,
                                    sync_workflow_runner=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.runner_dir = self.generate_tmp_runner_dir()
        self.runner_paths = list(glob.glob(
            '%s/**/*_wf_runner.py' % self.runner_dir, recursive=True))

    def generate_tmp_runner_dir(self):
        tmp_runner_dir = tempfile.mkdtemp(prefix='sync_wfr.test.')
        for i in range(3):
            runner_path = os.path.join(tmp_runner_dir,
                                       'runner_%s_wf_runner.py' % i)
            with open(runner_path, 'w') as f: f.write("some content")
        return tmp_runner_dir

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()
        super().tearDown()

    def test_falls_back_to_default_runner_dirs(self):
        utils.sync_workflow_runners()
        self.assertEqual(
            self.mocks['utils']['get_default_runner_dirs'].call_count, 1)

    def test_registers_new_runners(self):
        result = utils.sync_workflow_runners(runner_dirs=[self.runner_dir])
        mock_register_fn = self.mocks['utils']['add_workflow_runner']
        call_args_list = mock_register_fn.call_args_list
        expected_call_args_list = [call(runner_path=runner_path)
                                   for runner_path in self.runner_paths]
        self.assertEqual(
            sorted(call_args_list, key=lambda c: c[1]['runner_path']),
            sorted(expected_call_args_list, key=lambda c: c[2]['runner_path']),
        )
        self.assertEqual(result['created'],
                         {path: mock_register_fn.return_value
                          for path in self.runner_paths})

    def test_handles_registration_errors(self):
        mock_register_fn = self.mocks['utils']['add_workflow_runner']
        mock_register_fn.side_effect = Exception()
        result = utils.sync_workflow_runners(runner_dirs=[self.runner_dir])
        self.assertEqual(
            result['errors'],
            {path: mock_register_fn.side_effect for path in self.runner_paths}
        )

    def test_marks_missing_runners(self):
        missing_runner = WorkflowRunner.objects.create(path='/a/missing/path')
        result = utils.sync_workflow_runners(runner_dirs=[self.runner_dir])
        self.assertEqual(
            set(WorkflowRunner.objects.exclude(error__isnull=True)),
            set([missing_runner])
        )
        self.assertTrue(isinstance(result['errors'][missing_runner.path],
                                   utils.MissingRunnerFileError))

    def test_syncs_existing_runners(self):
        existing_runners = [WorkflowRunner.objects.create(path=path)
                            for path in self.runner_paths]
        result = utils.sync_workflow_runners(runner_dirs=[self.runner_dir])
        mock_sync_fn = self.mocks['utils']['sync_workflow_runner']
        call_args_list = mock_sync_fn.call_args_list
        expected_call_args_list = [call(runner_path=runner.path)
                                   for runner in existing_runners]
        self.assertEqual(
            sorted(call_args_list, key=lambda c: c[1]['runner_path']),
            sorted(expected_call_args_list, key=lambda c: c[2]['runner_path']),
        )
        self.assertEqual(
            result['synced'],
            {runner.path: mock_sync_fn.return_value
             for runner in existing_runners}
        )

    def test_handles_sync_errors(self):
        existing_runners = [WorkflowRunner.objects.create(path=path)
                            for path in self.runner_paths]
        mock_sync_fn = self.mocks['utils']['sync_workflow_runner']
        mock_sync_fn.side_effect = Exception()
        result = utils.sync_workflow_runners(runner_dirs=[self.runner_dir])
        self.assertEqual(
            set(WorkflowRunner.objects.exclude(error__isnull=True)),
            set(existing_runners)
        )
        self.assertEqual(
            result['errors'],
            {runner.path: mock_sync_fn.side_effect
             for runner in existing_runners}
        )
