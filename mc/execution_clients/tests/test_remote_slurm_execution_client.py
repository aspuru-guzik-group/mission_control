from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from ..remote_slurm_execution_client import RemoteSlurmExecutionClient


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.ssh_client = MagicMock()
        self.rsclient = RemoteSlurmExecutionClient(
            ssh_client=self.ssh_client
        )
        self.slurm_execution_client = MagicMock()
        self.rsclient.slurm_execution_client = self.slurm_execution_client
        self.submission = defaultdict(MagicMock)
        self.execution_meta = defaultdict(MagicMock)

class StartExecutionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.rsclient.upload_submission = MagicMock()

    def test_uploads_submission(self):
        self.rsclient.start_execution(submission=self.submission)
        self.assertEqual(self.rsclient.upload_submission.call_args,
                         call(submission=self.submission))

    def test_calls_slurm_execution_client_method_w_remote_dir(self):
        wrapped_fn = self.rsclient.slurm_execution_client.start_execution
        self.rsclient.start_execution(submission=self.submission)
        expected_modified_submission = {
            **self.submission, 
            'dir': self.rsclient.upload_submission.return_value['dir']
        }
        self.assertEqual(wrapped_fn.call_args,
                         call(submission=expected_modified_submission))

class UploadSubmissionDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.rsclient.generate_remote_dir_path = MagicMock()
        self.expected_remote_dir_path = \
                self.rsclient.generate_remote_dir_path.return_value
    
    def test_ensures_remote_parent_dir(self):
        self.rsclient.upload_submission(submission=self.submission)
        self.assertEqual(
            self.rsclient.ssh_client.run_process.call_args,
            call(cmd=['mkdir', '-p', self.rsclient.remote_workdir], check=True)
        )

    def test_calls_rsync(self):
        self.rsclient.upload_submission(submission=self.submission)
        self.assertEqual(
            self.rsclient.ssh_client.rsync_to_remote.call_args,
            call(local_src_path=self.submission['dir'] + '/',
                 remote_dest_path=self.expected_remote_dir_path,
                 flags=self.rsclient.rsync_flags)
        )

    def test_returns_remote_dir_meta(self):
        remote_dir_meta = self.rsclient.upload_submission(
            submission=self.submission)
        self.assertEqual(remote_dir_meta['dir'], self.expected_remote_dir_path)

class GetExecutionStateTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.rsclient.on_execution_completed = MagicMock()
        self.wrapped_get_execution_state = \
                self.rsclient.slurm_execution_client.get_execution_state

    def test_calls_slurm_execution_client_method(self):
        execution_state = self.rsclient.get_execution_state(
            execution_meta=self.execution_meta)
        self.assertEqual(self.wrapped_get_execution_state.call_args,
                         call(execution_meta=self.execution_meta))
        self.assertEqual(execution_state,
                         self.wrapped_get_execution_state.return_value)

    def test_calls_on_execution_completed_if_submission_not_running(self):
        self.wrapped_get_execution_state.return_value = {'run_status': ''}
        self.rsclient.get_execution_state(execution_meta=self.execution_meta)
        self.assertEqual(self.rsclient.on_execution_completed.call_count, 1)

    def test_does_not_call_on_execution_completed_if_submission_running(self):
        self.wrapped_get_execution_state.return_value = \
                {'run_status': 'RUNNING'}
        self.rsclient.get_execution_state(execution_meta=self.execution_meta)
        self.assertEqual(self.rsclient.on_execution_completed.call_count, 0)

class OnExecutionCompletedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.rsclient.download_completed_dir = MagicMock()
        self.rsclient.on_execution_completed(execution_meta=self.execution_meta)

    def test_downloads_submission_dir(self):
        self.assertEqual(self.rsclient.download_completed_dir.call_args,
                         call(execution_meta=self.execution_meta))

class DownloadCompletedDirTestCase(BaseTestCase):
    def test_calls_rsync(self):
        self.rsclient.download_completed_dir(execution_meta=self.execution_meta)
        expected_remote_src_path = \
                self.execution_meta['remote_dir']['dir'] + '/'
        expected_local_dest_path = self.execution_meta['submission']['dir']
        self.assertEqual(
            self.rsclient.ssh_client.rsync_from_remote.call_args,
            call(remote_src_path=expected_remote_src_path,
                 local_dest_path=expected_local_dest_path,
                 flags=self.rsclient.rsync_flags)
        )

if __name__ == '__main__':
    unittest.main()
