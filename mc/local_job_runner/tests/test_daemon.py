import operator
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch

from ..daemon import LocalJobRunnerDaemon

class DaemonBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_job_client = MagicMock()
        self.mock_job_dir_factory = MagicMock()
        self.daemon = LocalJobRunnerDaemon(
            job_client=self.mock_job_client,
            job_dir_factory=self.mock_job_dir_factory)

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
        self.patcher = patch.multiple(
            self.daemon, fetch_candidate_job_specs=DEFAULT,
            process_candidate_job_spec=DEFAULT, process_running_jobs=DEFAULT,
            process_transferring_jobs=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_fetches_candidate_job_spec(self):
        self.daemon.tick()
        self.assertTrue(self.mocks['fetch_candidate_job_specs'].called)

    def test_processes_candidate_job_specs(self):
        self.daemon.max_running_jobs = 3
        extra_jobs = 2
        job_specs = [{'uuid': i}
                     for i in range(self.daemon.max_running_jobs + extra_jobs)]
        self.mocks['fetch_candidate_job_specs'].return_value = job_specs
        self.daemon.tick()
        expected_call_args = [call(job_spec=job_spec)
                              for job_spec in job_specs[:-1 * extra_jobs]]
        self.assertTrue(self.mocks['process_candidate_job_spec'].call_args_list,
                        expected_call_args)

    def test_processes_running_jobs(self):
        self.daemon.tick()
        self.assertEqual(self.mocks['process_running_jobs'].call_count, 1)

    def test_processes_transferring_jobs(self):
        self.daemon.tick()
        self.assertEqual(self.mocks['process_transferring_jobs'].call_count, 1)

class FetchCandidateJobsTestCase(DaemonBaseTestCase):
    def test_fetch_candidate_job_specs(self):
        self.daemon.fetch_candidate_job_specs()
        self.assertEqual(self.mock_job_client.fetch_jobs.call_count, 1)

class ProcessCandidateJobSpecTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon, build_job_dir=DEFAULT,
                                      run_job=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_claimable_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        self.mock_job_client.claim_jobs.return_value = {job_spec['uuid']: True}
        self.daemon.process_candidate_job_spec(job_spec=job_spec)
        expected_job_dir_meta = self.mocks['build_job_dir'].return_value
        expected_job_proc_meta = self.mocks['run_job'].return_value
        self.assertEqual(self.mocks['build_job_dir'].call_count, 1)
        expected_partial_job = {
            'job_spec': job_spec,
            'dir': expected_job_dir_meta,
        }
        self.assertEqual(self.mocks['run_job'].call_args,
                         call(job=expected_partial_job))
        expected_full_job = {**expected_partial_job,
                             'proc': expected_job_proc_meta}
        self.assertEqual(self.daemon.running_jobs[job_spec['uuid']],
                         expected_full_job)

    def test_unclaimable_job(self):
        mock_job = {'uuid': 'abcd'}
        self.mock_job_client.claim_jobs.return_value = {mock_job['uuid']: False}
        self.daemon.process_candidate_job_spec(mock_job)
        self.assertEqual(self.mocks['run_job'].call_count, 0)

class BuildJobDirTestCase(DaemonBaseTestCase):
    def test_build_job_dir(self):
        job_spec = {'uuid': 'abcd'}
        self.daemon.build_job_dir(job_spec=job_spec)
        self.assertTrue(self.daemon.job_dir_factory.build_job_dir.call_args,
                        call(job_spec=job_spec))

class RunJobTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch('subprocess.Popen')
        self.mocks = {'subprocess.Popen': self.patcher.start()}

    def tearDown(self):
        self.patcher.stop()

    def test_run_job(self):
        # calls Popen w/ cd job dir, command.
        job = {
            'job_spec': {'uuid': 'abcd'},
            'dir': {'dir': 'mock_dir', 'entrypoint': 'job.sh'},
        }
        proc_meta = self.daemon.run_job(job=job)
        expected_cmd = 'cd {dir}; {entrypoint}'.format(**job['dir'])
        expected_popen_call = call(expected_cmd, shell=True)
        self.assertEqual(self.mocks['subprocess.Popen'].call_args,
                         expected_popen_call)
        expected_proc_meta = {
            'pid': self.mocks['subprocess.Popen'].return_value.pid}
        self.assertEqual(proc_meta, expected_proc_meta)

class ProcessRunningJobsTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.completed_jobs = [{'job_spec': {'uuid': 'completed_%s' % i}}
                               for i in range(3)]
        self.incompleted_jobs = [{'job_spec': {'uuid': 'incompleted_%s' % i}}
                                 for i in range(3)]
        self.daemon.running_jobs = {
            j['job_spec']['uuid']: j
            for j in self.completed_jobs + self.incompleted_jobs
        }
        self.mock_job_states = {
            **{j['job_spec']['uuid']: None
               for j in self.completed_jobs},
            **{j['job_spec']['uuid']: {'mock': 'state'}
               for j in self.incompleted_jobs}
        }
        self.patcher = patch.multiple(self.daemon,
                                      get_running_job_states=DEFAULT,
                                      process_completed_job=DEFAULT)
        self.mocks = self.patcher.start()
        self.mocks['get_running_job_states'].return_value = self.mock_job_states

    def tearDown(self):
        self.patcher.stop()

    def test_processes_completed_jobs(self):
        self.daemon.process_running_jobs()
        expected_calls = [call(job=job) for job in self.completed_jobs]
        calls = self.mocks['process_completed_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['job_spec']['uuid']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)

class ProcessCompletedJobTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon, transfer_job=DEFAULT)
        self.mocks = self.patcher.start()
        self.job_uuid = 'abcd'
        self.completed_job = {'job_spec': {'uuid': self.job_uuid}}
        self.daemon.running_jobs[self.job_uuid] = self.completed_job

    def tearDown(self):
        self.patcher.stop()

    def test_transfers_completed_jobs(self):
        self.daemon.process_completed_job(job=self.completed_job)
        self.assertEqual(self.mocks['transfer_job'].call_args_list,
                         [call(job=self.completed_job)])
        self.assertTrue(self.job_uuid not in self.daemon.running_jobs)
        self.assertEqual(self.daemon.transferring_jobs[self.job_uuid],
                         {**self.completed_job,
                          'transfer': self.mocks['transfer_job'].return_value})

class GetRunningJobStatesTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch('psutil.Process')
        self.mocks = {'psutil.Process': self.patcher.start()}

        self.daemon.running_jobs = {}
        for i in range(3):
            job_uuid = 'job_uuid_%s' % i
            self.daemon.running_jobs[job_uuid] = {
                'job': {'uuid': job_uuid},
                'dir': {'dir': 'dir_%s' % job_uuid},
                'proc': {'pid': i},
            }

    def tearDown(self):
        self.patcher.stop()

    def test_get_running_job_states(self):
        job_states = self.daemon.get_running_job_states()
        expected_psutil_calls = [call(pid=j['proc']['pid'])
                                 for j in self.daemon.running_jobs.values()]
        self.assertEqual(self.mocks['psutil.Process'].call_args_list,
                         expected_psutil_calls)
        mock_job_state = operator.attrgetter(
            'return_value.as_dict.return_value')(self.mocks['psutil.Process'])
        expected_job_states = {job_uuid: mock_job_state
                               for job_uuid in self.daemon.running_jobs.keys()}
        self.assertEqual(job_states, expected_job_states)

    def test_missing_process(self):
        self.fail("ADD TEST LOGIC HERE FOR WHEN psutil.Process fails")

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
