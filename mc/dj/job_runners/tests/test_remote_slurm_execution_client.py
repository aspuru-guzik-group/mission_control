import os
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from ..remote_slurm_execution_client import RemoteSlurmExecutionClient


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.ssh_client = MagicMock()
        self.rsclient = RemoteSlurmExecutionClient(
            ssh_client=self.ssh_client
        )
        self.slurm_execution_client = MagicMock()
        self.rsclient.slurm_execution_client = self.slurm_execution_client
        self.job = MagicMock()

    def tearDown(self):
        self.stop_patchers()

    def start_patchers(self):
        if hasattr(self, 'patchers'):
            self.mocks = {key: patcher.start()
                          for key, patcher in self.patchers.items()}

    def stop_patchers(self):
        if hasattr(self, 'patchers'):
            for patcher in self.patchers.values(): patcher.stop()

class StartExecutionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'rsclient': patch.multiple(self.rsclient, upload_job=DEFAULT)
        }
        self.start_patchers()

    def test_uploads_job(self):
        self.rsclient.start_execution(job=self.job)
        self.assertEqual(self.rsclient.upload_job.call_args,
                         call(job=self.job))

    def test_calls_slurm_execution_client_method_w_remote_dir(self):
        wrapped_fn = self.rsclient.slurm_execution_client.start_execution
        self.rsclient.start_execution(job=self.job)
        expected_job_arg = {
            **self.job, 
            'submission': {
                **self.job['submission'],
                'dir': self.rsclient.upload_job.return_value['dir']
            }
        }
        self.assertEqual(wrapped_fn.call_args, call(job=expected_job_arg))

class UploadJobDirTestCase(BaseTestCase):
    def test_ensures_remote_parent_dir(self):
        self.rsclient.upload_job(job=self.job)
        self.assertEqual(
            self.rsclient.ssh_client.run_process.call_args,
            call(cmd=['mkdir', '-p', self.rsclient.remote_workdir], check=True)
        )

    def test_calls_rsync(self):
        self.rsclient.upload_job(job=self.job)
        self.assertEqual(
            self.rsclient.ssh_client.rsync_to_remote.call_args,
            call(local_src_path=self.job['submission']['dir'] + '/',
                 remote_dest_path=os.path.join(self.rsclient.remote_workdir,
                                               self.job['uuid']),
                 flags='-a')
        )

    def test_returns_remote_dir_meta(self):
        remote_dir_meta = self.rsclient.upload_job(job=self.job)
        self.assertEqual(
            remote_dir_meta['dir'],
            os.path.join(self.rsclient.remote_workdir, self.job['uuid'])
        )

class GetExecutionStateTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'rsclient': patch.multiple(self.rsclient, on_job_completed=DEFAULT)
        }
        self.wrapped_fn = self.rsclient.slurm_execution_client\
                .get_execution_state
        self.start_patchers()

    def test_calls_slurm_execution_client_method(self):
        execution_state = self.rsclient.get_execution_state(job=self.job)
        self.assertEqual(self.wrapped_fn.call_args, call(job=self.job))
        self.assertEqual(execution_state, self.wrapped_fn.return_value)

    def test_calls_on_completed_if_job_not_executing(self):
        self.wrapped_fn.return_value = {'executing': False}
        self.rsclient.get_execution_state(job=self.job)
        self.assertEqual(self.rsclient.on_job_completed.call_count, 1)

    def test_does_not_call_on_completed_if_job_executing(self):
        self.wrapped_fn.return_value = {'executing': True}
        self.rsclient.get_execution_state(job=self.job)
        self.assertEqual(self.rsclient.on_job_completed.call_count, 0)

class OnCompletedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'rsclient': patch.multiple(self.rsclient, download_job=DEFAULT),
        }
        self.start_patchers()
        self.rsclient.on_job_completed(job=self.job)

    def test_downloads_job_dir(self):
        self.assertEqual(self.rsclient.download_job.call_args,
                         call(job=self.job))

    def test_sets_completed_dir_in_execution_meta(self):
        self.assertEqual(self.job['execution']['completed_dir'],
                         self.job['dir']['dir'])

class DownloadJobTestCase(BaseTestCase):
    def test_calls_rsync(self):
        self.rsclient.download_job(job=self.job)
        expected_remote_src_path = (
            self.job['execution']['remote_dir']['dir'] + '/')
        self.assertEqual(
            self.rsclient.ssh_client.rsync_from_remote.call_args,
            call(remote_src_path=expected_remote_src_path,
                 local_dest_path=self.job['dir']['dir'], flags='-a')
        )

if __name__ == '__main__':
    unittest.main()
