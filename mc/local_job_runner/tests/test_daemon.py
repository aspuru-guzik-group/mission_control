import operator
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch

from ..daemon import LocalJobRunnerDaemon

class DaemonBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job_client = MagicMock()
        self.job_dir_factory = MagicMock()
        self.transfer_client = MagicMock()
        self.daemon = LocalJobRunnerDaemon(
            job_client=self.job_client,
            job_dir_factory=self.job_dir_factory,
            transfer_client=self.transfer_client)

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
            process_candidate_job_spec=DEFAULT, process_executing_jobs=DEFAULT,
            process_transferring_jobs=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_fetches_candidate_job_spec(self):
        self.daemon.tick()
        self.assertTrue(self.mocks['fetch_candidate_job_specs'].called)

    def test_processes_candidate_job_specs(self):
        self.daemon.max_executing_jobs = 3
        extra_jobs = 2
        job_specs = [
            {'uuid': i}
            for i in range(self.daemon.max_executing_jobs + extra_jobs)
        ]
        self.mocks['fetch_candidate_job_specs'].return_value = job_specs
        self.daemon.tick()
        expected_call_args = [call(job_spec=job_spec)
                              for job_spec in job_specs[:-1 * extra_jobs]]
        self.assertTrue(self.mocks['process_candidate_job_spec'].call_args_list,
                        expected_call_args)

    def test_processes_executing_jobs(self):
        self.daemon.tick()
        self.assertEqual(self.mocks['process_executing_jobs'].call_count, 1)

    def test_processes_transferring_jobs(self):
        self.daemon.tick()
        self.assertEqual(self.mocks['process_transferring_jobs'].call_count, 1)

class FetchCandidateJobsTestCase(DaemonBaseTestCase):
    def test_fetch_candidate_job_specs(self):
        self.daemon.fetch_candidate_job_specs()
        self.assertEqual(self.job_client.fetch_jobs.call_count, 1)

class ProcessCandidateJobSpecTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon, build_job_dir=DEFAULT,
                                      start_job_execution=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_claimable_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {job_spec['uuid']: True}
        self.daemon.process_candidate_job_spec(job_spec=job_spec)
        expected_job_dir_meta = self.mocks['build_job_dir'].return_value
        expected_job_proc_meta = self.mocks['start_job_execution'].return_value
        self.assertEqual(self.mocks['build_job_dir'].call_count, 1)
        expected_partial_job = {
            'key': job_spec['uuid'],
            'job_spec': job_spec,
            'dir': expected_job_dir_meta,
        }
        self.assertEqual(self.mocks['start_job_execution'].call_args,
                         call(job=expected_partial_job))
        expected_full_job = {**expected_partial_job,
                             'proc': expected_job_proc_meta}
        self.assertEqual(
            self.daemon.executing_jobs[expected_partial_job['key']],
            expected_full_job)

    def test_unclaimable_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {job_spec['uuid']: False}
        self.daemon.process_candidate_job_spec(job_spec)
        self.assertEqual(self.mocks['start_job_execution'].call_count, 0)

class BuildJobDirTestCase(DaemonBaseTestCase):
    def test_build_job_dir(self):
        job_spec = {'uuid': 'abcd'}
        self.daemon.build_job_dir(job_spec=job_spec)
        self.assertTrue(self.daemon.job_dir_factory.build_job_dir.call_args,
                        call(job_spec=job_spec))

class StartJobExecutionTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch('subprocess.Popen')
        self.mocks = {'subprocess.Popen': self.patcher.start()}

    def tearDown(self):
        self.patcher.stop()

    def test_start_job_execution(self):
        # calls Popen w/ cd job dir, command.
        job = {
            'key': 'abcd',
            'dir': {'dir': 'mock_dir', 'entrypoint': 'job.sh'},
        }
        proc_meta = self.daemon.start_job_execution(job=job)
        expected_cmd = 'cd {dir}; {entrypoint}'.format(**job['dir'])
        expected_popen_call = call(expected_cmd, shell=True)
        self.assertEqual(self.mocks['subprocess.Popen'].call_args,
                         expected_popen_call)
        expected_proc_meta = {
            'pid': self.mocks['subprocess.Popen'].return_value.pid}
        self.assertEqual(proc_meta, expected_proc_meta)

class ProcessExecutingJobsTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.executed_jobs = [{'key': 'executed_%s' % i} for i in range(3)]
        self.executing_jobs = [{'key': 'executing_%s' % i} for i in range(3)]
        self.daemon.executing_jobs = {
            j['key']: j for j in self.executed_jobs + self.executing_jobs
        }
        self.mock_job_execution_states = {
            **{j['key']: None for j in self.executed_jobs},
            **{j['key']: {'mock': 'state'} for j in self.executing_jobs}
        }
        self.patcher = patch.multiple(self.daemon,
                                      get_job_execution_states=DEFAULT,
                                      process_executed_job=DEFAULT)
        self.mocks = self.patcher.start()
        self.mocks['get_job_execution_states'].return_value = \
                self.mock_job_execution_states

    def tearDown(self):
        self.patcher.stop()

    def test_process_executing_jobs(self):
        self.daemon.process_executing_jobs()
        expected_calls = [call(job=job) for job in self.executed_jobs]
        calls = self.mocks['process_executed_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['key']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)

class ProcessCompletedJobTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon, start_job_transfer=DEFAULT)
        self.mocks = self.patcher.start()
        self.job_key = 'abcd'
        self.executed_job = {'key': self.job_key}
        self.daemon.executing_jobs[self.job_key] = self.executed_job

    def tearDown(self):
        self.patcher.stop()

    def test_transfers_executed_jobs(self):
        self.daemon.process_executed_job(job=self.executed_job)
        self.assertEqual(self.mocks['start_job_transfer'].call_args_list,
                         [call(job=self.executed_job)])
        self.assertTrue(self.job_key not in self.daemon.executing_jobs)
        self.assertEqual(
            self.daemon.transferring_jobs[self.job_key],
            {**self.executed_job,
             'transfer': self.mocks['start_job_transfer'].return_value})

class GetJobRunStatesTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch('psutil.Process')
        self.mocks = {'psutil.Process': self.patcher.start()}

        self.daemon.executing_jobs = {}
        for i in range(3):
            job_key = 'job_key_%s' % i
            self.daemon.executing_jobs[job_key] = {'key': job_key,
                                                 'proc': {'pid': i}}

    def tearDown(self):
        self.patcher.stop()

    def test_get_job_execution_states(self):
        job_execution_states = self.daemon.get_job_execution_states()
        expected_psutil_calls = [call(pid=j['proc']['pid'])
                                 for j in self.daemon.executing_jobs.values()]
        self.assertEqual(self.mocks['psutil.Process'].call_args_list,
                         expected_psutil_calls)
        mock_job_execution_state = operator.attrgetter(
            'return_value.as_dict.return_value')(self.mocks['psutil.Process'])
        expected_job_states = {job_key: mock_job_execution_state
                               for job_key in self.daemon.executing_jobs.keys()}
        self.assertEqual(job_execution_states, expected_job_states)

class TransferJobTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_key = 'abcd'
        self.job = {'key': self.job_key}

    def test_start_transfer_job(self):
        self.daemon.start_job_transfer(job=self.job)
        self.assertEqual(self.transfer_client.start_transfer.call_args,
                         call(job=self.job))


class ProcessTransferringJobsTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.transferred_jobs = [{'key': 'transferred_%s' % i}
                                 for i in range(3)]
        self.transferring_jobs = [{'key': 'transferring_%s' % i}
                                  for i in range(3)]
        self.daemon.transferring_jobs = {
            j['key']: j for j in self.transferred_jobs + self.transferring_jobs
        }
        self.mock_job_transfer_states = {
            **{j['key']: None for j in self.transferred_jobs},
            **{j['key']: {'mock': 'state'} for j in self.transferring_jobs}
        }
        self.patcher = patch.multiple(self.daemon,
                                      get_job_transfer_states=DEFAULT,
                                      process_transferred_job=DEFAULT)
        self.mocks = self.patcher.start()
        self.mocks['get_job_transfer_states'].return_value = \
                self.mock_job_transfer_states

    def tearDown(self):
        self.patcher.stop()

    def test_process_transferring_jobs(self):
        self.daemon.process_transferring_jobs()
        expected_calls = [call(job=job) for job in self.transferred_jobs]
        calls = self.mocks['process_transferred_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['key']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)


class GetJobTransferStatesTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.daemon.transferring_jobs = {}
        for i in range(3):
            job_key = 'job_%s' % i
            self.daemon.transferring_jobs[job_key] = {
                'key': job_key,
                'transfer': {'key': 'transfer_%s' % job_key},
            }

    def test_get_job_transfer_states(self):
        job_transfer_states = self.daemon.get_job_transfer_states()
        expected_calls = [call(key=j['transfer']['key'])
                          for j in self.daemon.transferring_jobs.values()]
        self.assertEqual(self.transfer_client.get_transfer_state.call_args_list,
                         expected_calls)
        expected_job_states = {
            job_key: self.transfer_client.get_transfer_state.return_value
            for job_key in self.daemon.transferring_jobs.keys()
        }
        self.assertEqual(job_transfer_states, expected_job_states)

class ProcessTransferredJobTestCase(DaemonBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.daemon, update_job_spec=DEFAULT)
        self.mocks = self.patcher.start()
        self.job_key = 'abcd'
        self.transferred_job = {
            'key': self.job_key,
            'job_spec': {'uuid': 'uuid_%s' % self.job_key},
            'transfer': {'junk': 'junk'}
        }
        self.daemon.transferring_jobs[self.job_key] = self.transferred_job

    def tearDown(self):
        self.patcher.stop()

    def test_process_transferred_job(self):
        self.daemon.process_transferred_job(job=self.transferred_job)
        self.assertEqual(
            self.mocks['update_job_spec'].call_args,
            call(job_spec=self.transferred_job['job_spec'], updates={
                'status': self.job_client.statuses.TRANSFERRED,
                'transfer_meta': self.transferred_job['transfer']
            }))
        self.assertTrue(self.job_key not in self.daemon.transferring_jobs)

class UpdateJobSpecTestCase(DaemonBaseTestCase):
    def test_update_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        updates = {'pie': 'blueberry', 'meat': 'beef'}
        self.daemon.update_job_spec(job_spec=job_spec, updates=updates)
        self.assertEqual(self.job_client.update_job.call_args, 
                         call(uuid=job_spec['uuid'], updates=updates))

if __name__ == '__main__':
    unittest.main()
