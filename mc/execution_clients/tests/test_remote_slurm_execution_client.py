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
    def test_wraps_slurm_execution_client_method(self):
        wrapped_get_execution_state = \
                self.rsclient.slurm_execution_client.get_execution_state
        execution_state = self.rsclient.get_execution_state(
            execution_meta=self.execution_meta)
        self.assertEqual(wrapped_get_execution_state.call_args,
                         call(execution_meta=self.execution_meta))
        self.assertEqual(execution_state,
                         wrapped_get_execution_state.return_value)

if __name__ == '__main__':
    unittest.main()
