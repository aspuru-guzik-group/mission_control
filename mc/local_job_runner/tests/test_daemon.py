import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch

from ..daemon import LocalJobRunnerDaemon

class DaemonBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_job_client = MagicMock()
        self.daemon = LocalJobRunnerDaemon(job_client=self.mock_job_client)

class RunTestCase(DaemonBaseTestCase):
    def test_run(self):
        intervals = [1, 2]
        self.daemon.interval = 10
        with patch('time.sleep') as mock_sleep:
            for interval in intervals:
                self.daemon.interval = interval
                self.daemon.run(ntimes=1)
            self.assertEqual(mock_sleep.call_args_list,
                             [call(interval) for interval in intervals])

class TickTestCase(DaemonBaseTestCase):

    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon,
                                      fetch_candidate_jobs=DEFAULT,
                                      process_job=DEFAULT,
                                      poll_running_jobs=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_fetches_jobs(self):
        self.daemon.tick()
        self.assertTrue(self.mocks['fetch_candidate_jobs'].called)

    def test_processes_jobs(self):
        self.daemon.max_concurrent_jobs = 3
        extra_jobs = 2
        mock_jobs = [i for i in range(self.daemon.max_concurrent_jobs +
                                      extra_jobs)]
        self.mocks['fetch_candidate_jobs'].return_value = mock_jobs
        self.daemon.tick()
        self.assertTrue(self.mocks['process_job'].call_args_list,
                        [call(mock_job)
                         for mock_job in mock_jobs[:-1 * extra_jobs]])

    def test_polls_running_jobs(self):
        self.daemon.tick()
        self.assertEqual(self.mocks['poll_running_jobs'].call_count, 1)

class FetchCandidateJobsTestCase(DaemonBaseTestCase):
    def test_fetch_candidate_jobs(self):
        self.daemon.fetch_candidate_jobs()
        self.assertEqual(self.mock_job_client.fetch_jobs.call_count, 1)

class ProcessJobTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon,
                                      build_job_dir=DEFAULT,
                                      run_job=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_claimable_job(self):
        mock_job = {'uuid': 'abcd'}
        self.mock_job_client.claim_jobs.return_value = {mock_job['uuid']: True}
        self.daemon.process_job(mock_job)
        expected_job_dir = self.mocks['build_job_dir'].return_value
        expected_job_proc_meta = self.mocks['run_job'].return_value
        self.assertEqual(self.mocks['build_job_dir'].call_count, 1)
        self.assertEqual(self.mocks['run_job'].call_args,
                         call(job=mock_job, job_dir=expected_job_dir))
        self.assertEqual(self.daemon.running_jobs[mock_job['uuid']],
                         {'job': mock_job, 'run_dir': expected_job_dir,
                          'proc': expected_job_proc_meta})

    def test_unclaimable_job(self):
        mock_job = {'uuid': 'abcd'}
        self.mock_job_client.claim_jobs.return_value = {mock_job['uuid']: False}
        self.daemon.process_job(mock_job)
        self.assertEqual(self.mocks['run_job'].call_count, 0)

class RunJobTestCase(DaemonBaseTestCase):
    def test_run_job(self):
        self.fail()

class BuildJobDirTestCase(DaemonBaseTestCase):
    def test_build_job_dir(self):
        # Builds job dir.
        self.fail()

class PollRunningJobsTestCase(DaemonBaseTestCase):
    def test_poll_jobs(self):
        pass

class UploadCompletedJobTestCase(DaemonBaseTestCase):
    def test_upload_completed_job(self):
        # Posts job .tgz to endpoint.
        self.fail()

class UpdateJobStatusTestCase(DaemonBaseTestCase):
    def test_update_job_status(self):
        # Updates job status.
        self.fail()

if __name__ == '__main__':
    unittest.main()
