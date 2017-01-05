import logging
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from ..base_job_runner import BaseJobRunner


class JobRunnerBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job_spec_client = MagicMock()
        self.job_dir_factory = MagicMock()
        self.execution_client = MagicMock()
        self.transfer_client = MagicMock()
        self.runner = BaseJobRunner(
            execution_client=self.execution_client,
            job_spec_client=self.job_spec_client,
            job_dir_factory=self.job_dir_factory,
            transfer_client=self.transfer_client)

class RunTestCase(JobRunnerBaseTestCase):
    def test_run(self):
        tick_intervals = [1, 2]
        self.runner.tick_interval = 10
        with patch('time.sleep') as mock_sleep:
            for tick_interval in tick_intervals:
                self.runner.tick_interval = tick_interval
                self.runner.run(ntimes=1)
            self.assertEqual(
                mock_sleep.call_args_list,
                [call(tick_interval) for tick_interval in tick_intervals])

class TickTestCase(JobRunnerBaseTestCase):

    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(
            self.runner, fetch_claimable_job_specs=DEFAULT,
            process_claimable_job_spec=DEFAULT, process_executing_jobs=DEFAULT,
            process_transferring_jobs=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_fetches_claimable_job_spec(self):
        self.runner.tick()
        self.assertTrue(self.mocks['fetch_claimable_job_specs'].called)

    def test_processes_claimable_job_specs(self):
        self.runner.max_executing_jobs = 3
        extra_jobs = 2
        job_specs = [
            {'uuid': i}
            for i in range(self.runner.max_executing_jobs + extra_jobs)
        ]
        self.mocks['fetch_claimable_job_specs'].return_value = job_specs
        self.runner.tick()
        expected_call_args = [call(job_spec=job_spec)
                              for job_spec in job_specs[:-1 * extra_jobs]]
        self.assertTrue(self.mocks['process_claimable_job_spec'].call_args_list,
                        expected_call_args)

    def test_processes_executing_jobs(self):
        self.runner.tick()
        self.assertEqual(self.mocks['process_executing_jobs'].call_count, 1)

    def test_processes_transferring_jobs(self):
        self.runner.tick()
        self.assertEqual(self.mocks['process_transferring_jobs'].call_count, 1)

class FetchClaimableJobsTestCase(JobRunnerBaseTestCase):
    def test_fetch_claimable_job_specs(self):
        self.runner.fetch_claimable_job_specs()
        self.assertEqual(
            self.job_spec_client.fetch_claimable_job_specs.call_count, 1)

class ProcessClaimableJobSpecTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.runner, build_job_dir=DEFAULT,
                                      start_job_execution=DEFAULT,
                                      update_job_spec=DEFAULT)
        self.mocks = self.patcher.start()

    def tearDown(self):
        logging.disable(logging.NOTSET)
        self.patcher.stop()
        super().tearDown()

    def test_claimable_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        mock_claimed_spec = MagicMock()
        self.job_spec_client.claim_job_specs.return_value = {
            job_spec['uuid']: mock_claimed_spec}
        self.runner.process_claimable_job_spec(job_spec=job_spec)
        expected_job_dir_meta = self.mocks['build_job_dir'].return_value
        expected_execution_meta = self.mocks['start_job_execution'].return_value
        self.assertEqual(self.mocks['build_job_dir'].call_count, 1)
        expected_partial_job = {
            'key': mock_claimed_spec['uuid'],
            'job_spec': mock_claimed_spec,
            'dir': expected_job_dir_meta,
        }
        self.assertEqual(self.mocks['start_job_execution'].call_args,
                         call(job=expected_partial_job))
        expected_full_job = {**expected_partial_job,
                             'execution': expected_execution_meta}
        self.assertEqual(
            self.runner.executing_jobs[expected_partial_job['key']],
            expected_full_job)

    def test_unclaimable_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        self.job_spec_client.claim_job_specs.return_value = {
            job_spec['uuid']: None}
        self.runner.process_claimable_job_spec(job_spec)
        self.assertEqual(self.mocks['start_job_execution'].call_count, 0)

    def test_handles_start_execution_exception(self):
        logging.disable(logging.ERROR) # silently throw exception.
        job_spec = MagicMock()
        mock_claimed_spec = MagicMock()
        self.job_spec_client.claim_job_specs.return_value = {
            job_spec['uuid']: mock_claimed_spec}
        exception = Exception("some exception")
        self.mocks['start_job_execution'].side_effect = exception
        self.runner.process_claimable_job_spec(job_spec)
        self.assertEqual(
            self.mocks['update_job_spec'].call_args,
            call(job_spec=mock_claimed_spec, updates={'status': 'Failed'}))

class BuildJobDirTestCase(JobRunnerBaseTestCase):
    def test_build_job_dir(self):
        job_spec = {'uuid': 'abcd'}
        self.runner.build_job_dir(job_spec=job_spec)
        self.assertEqual(
            self.runner.job_dir_factory.build_dir_for_spec.call_args,
            call(job_spec=job_spec))

class StartJobExecutionTestCase(JobRunnerBaseTestCase):
    def test_start_job_execution(self):
        job = {'key': 'abcd'}
        execution_meta = self.runner.start_job_execution(job=job)
        self.assertEqual(self.execution_client.start_execution.call_args,
                         call(job=job))
        self.assertEqual(execution_meta,
                         self.execution_client.start_execution.return_value)

class ProcessExecutingJobsTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.executed_jobs = [{'key': 'executed_%s' % i} for i in range(3)]
        self.executing_jobs = [{'key': 'executing_%s' % i} for i in range(3)]
        self.runner.executing_jobs = {
            j['key']: j for j in self.executed_jobs + self.executing_jobs
        }
        self.mock_job_execution_states = {
            **{j['key']: None for j in self.executed_jobs},
            **{j['key']: {'executing': True} for j in self.executing_jobs}
        }
        self.patcher = patch.multiple(self.runner,
                                      get_job_execution_states=DEFAULT,
                                      process_executed_job=DEFAULT)
        self.mocks = self.patcher.start()
        self.mocks['get_job_execution_states'].return_value = \
                self.mock_job_execution_states

    def tearDown(self):
        self.patcher.stop()

    def test_process_executing_jobs(self):
        self.runner.process_executing_jobs()
        expected_calls = [call(job=job) for job in self.executed_jobs]
        calls = self.mocks['process_executed_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['key']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)

class GetJobExecutionStatesTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.executing_jobs = {}
        for i in range(3):
            job_key = 'job_key_%s' % i
            self.runner.executing_jobs[job_key] = {'key': job_key}

    def test_get_job_execution_states(self):
        job_execution_states = self.runner.get_job_execution_states()
        expected_get_job_execution_state_calls = [
            call(job=job) for job in self.runner.executing_jobs.values()
        ]
        self.assertEqual(
            self.execution_client.get_execution_state.call_args_list,
            expected_get_job_execution_state_calls)
        expected_execution_states = {
            job['key']: self.execution_client.get_execution_state.return_value
            for job in self.runner.executing_jobs.values()
        }
        self.assertEqual(job_execution_states, expected_execution_states)

class ProcessExecutedJobTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.runner, start_job_transfer=DEFAULT)
        self.mocks = self.patcher.start()
        self.job_key = 'abcd'
        self.executed_job = {'key': self.job_key}
        self.runner.executing_jobs[self.job_key] = self.executed_job

    def tearDown(self):
        self.patcher.stop()

    def test_transfers_executed_jobs(self):
        self.runner.process_executed_job(job=self.executed_job)
        self.assertEqual(self.mocks['start_job_transfer'].call_args_list,
                         [call(job=self.executed_job)])
        self.assertTrue(self.job_key not in self.runner.executing_jobs)
        self.assertEqual(
            self.runner.transferring_jobs[self.job_key],
            {**self.executed_job,
             'transfer': self.mocks['start_job_transfer'].return_value})

class TransferJobTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_key = 'abcd'
        self.job = {'key': self.job_key}

    def test_start_transfer_job(self):
        self.runner.start_job_transfer(job=self.job)
        self.assertEqual(self.transfer_client.start_transfer.call_args,
                         call(job=self.job))


class ProcessTransferringJobsTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.transferred_jobs = [{'key': 'transferred_%s' % i}
                                 for i in range(3)]
        self.transferring_jobs = [{'key': 'transferring_%s' % i}
                                  for i in range(3)]
        self.runner.transferring_jobs = {
            j['key']: j for j in self.transferred_jobs + self.transferring_jobs
        }
        self.mock_job_transfer_states = {
            **{j['key']: None for j in self.transferred_jobs},
            **{j['key']: {'transferring': True} for j in self.transferring_jobs}
        }
        self.patcher = patch.multiple(self.runner,
                                      get_job_transfer_states=DEFAULT,
                                      process_transferred_job=DEFAULT)
        self.mocks = self.patcher.start()
        self.mocks['get_job_transfer_states'].return_value = \
                self.mock_job_transfer_states

    def tearDown(self):
        self.patcher.stop()

    def test_process_transferring_jobs(self):
        self.runner.process_transferring_jobs()
        expected_calls = [call(job=job) for job in self.transferred_jobs]
        calls = self.mocks['process_transferred_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['key']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)


class GetJobTransferStatesTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.transferring_jobs = {}
        for i in range(3):
            job_key = 'job_%s' % i
            self.runner.transferring_jobs[job_key] = {
                'key': job_key,
                'transfer': {'key': 'transfer_%s' % job_key},
            }

    def test_get_job_transfer_states(self):
        job_transfer_states = self.runner.get_job_transfer_states()
        expected_calls = [call(job=j)
                          for j in self.runner.transferring_jobs.values()]
        self.assertEqual(self.transfer_client.get_transfer_state.call_args_list,
                         expected_calls)
        expected_job_states = {
            job_key: self.transfer_client.get_transfer_state.return_value
            for job_key in self.runner.transferring_jobs.keys()
        }
        self.assertEqual(job_transfer_states, expected_job_states)

class ProcessTransferredJobTestCase(JobRunnerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = patch.multiple(self.runner, update_job_spec=DEFAULT)
        self.mocks = self.patcher.start()
        self.job_key = 'abcd'
        self.transferred_job = {
            'key': self.job_key,
            'job_spec': {'uuid': 'uuid_%s' % self.job_key},
            'transfer': {'junk': 'junk'}
        }
        self.runner.transferring_jobs[self.job_key] = self.transferred_job

    def tearDown(self):
        self.patcher.stop()

    def test_process_transferred_job(self):
        self.runner.process_transferred_job(job=self.transferred_job)
        self.assertEqual(
            self.mocks['update_job_spec'].call_args,
            call(job_spec=self.transferred_job['job_spec'], updates={
                'status': self.job_spec_client.Statuses.Completed.name,
                'transfer_meta': self.transferred_job['transfer']
            }))
        self.assertTrue(self.job_key not in self.runner.transferring_jobs)

class UpdateJobSpecTestCase(JobRunnerBaseTestCase):
    def test_update_job_spec(self):
        job_spec = {'uuid': 'abcd'}
        updates = {'pie': 'blueberry', 'meat': 'beef'}
        self.runner.update_job_spec(job_spec=job_spec, updates=updates)
        self.assertEqual(self.job_spec_client.update_job_specs.call_args, 
                         call(updates_by_uuid={job_spec['uuid']: updates}))

if __name__ == '__main__':
    unittest.main()
