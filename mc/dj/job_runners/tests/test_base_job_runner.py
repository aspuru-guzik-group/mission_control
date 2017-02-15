from collections import defaultdict
import logging
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from ..base_job_runner import BaseJobRunner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.action_processor = MagicMock()
        self.execution_client = MagicMock()
        self.job_client = MagicMock()
        self.job_submission_factory = MagicMock()
        self.runner = BaseJobRunner(
            action_processor=self.action_processor,
            execution_client=self.execution_client,
            job_client=self.job_client,
            job_submission_factory=self.job_submission_factory)
        self.runner_methods_to_patch = []
        self.decorate_runner_methods_to_patch()
        if self.runner_methods_to_patch:
            self.patcher = patch.multiple(
                self.runner,
                **{method: DEFAULT for method in self.runner_methods_to_patch}
            )
            self.mocks = self.patcher.start()

    def decorate_runner_methods_to_patch(self): pass

    def tearDown(self):
        if hasattr(self, 'patcher'): self.patcher.stop()

class RunTestCase(BaseTestCase):
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

class TickTestCase(BaseTestCase):
    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['fetch_claimable_jobs',
                                             'process_claimable_job',
                                             'process_executing_jobs'])

    def test_fetches_claimable_job(self):
        self.runner.tick()
        self.assertTrue(self.mocks['fetch_claimable_jobs'].called)

    def test_processes_claimable_jobs(self):
        self.runner.max_executing_jobs = 3
        extra_jobs = 2
        jobs = [
            {'uuid': i}
            for i in range(self.runner.max_executing_jobs + extra_jobs)
        ]
        self.mocks['fetch_claimable_jobs'].return_value = jobs
        self.runner.tick()
        expected_call_args = [call(job=job)
                              for job in jobs[:-1 * extra_jobs]]
        self.assertTrue(self.mocks['process_claimable_job'].call_args_list,
                        expected_call_args)

    def test_processes_executing_jobs(self):
        self.runner.tick()
        self.assertEqual(self.mocks['process_executing_jobs'].call_count, 1)

class FetchClaimableJobsTestCase(BaseTestCase):
    def test_fetch_claimable_jobs(self):
        self.runner.fetch_claimable_jobs()
        self.assertEqual(
            self.job_client.fetch_claimable_jobs.call_count, 1)

class ProcessClaimableJobTestCase(BaseTestCase):
    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['build_job_submission',
                                             'start_job_execution',
                                             'update_job'])

    def tearDown(self):
        logging.disable(logging.NOTSET)
        super().tearDown()

    def test_claimable_job(self):
        job = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {job['uuid']: job}
        self.runner.process_claimable_job(job=job)
        expected_submission_meta = self.mocks['build_job_submission']\
                .return_value
        expected_execution_meta = self.mocks['start_job_execution'].return_value
        self.assertEqual(self.mocks['build_job_submission'].call_count, 1)
        self.assertEqual(self.mocks['start_job_execution'].call_args,
                         call(job=job))
        self.assertEqual(job['submission'], expected_submission_meta)
        self.assertEqual(job['execution'], expected_execution_meta)

    def test_unclaimable_job(self):
        job = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {
            job['uuid']: None}
        self.runner.process_claimable_job(job)
        self.assertEqual(self.mocks['start_job_execution'].call_count, 0)

    def test_handles_start_execution_exception(self):
        logging.disable(logging.ERROR) # silently throw exception.
        job = MagicMock()
        mock_claimed_job = MagicMock()
        self.job_client.claim_jobs.return_value = {
            job['uuid']: mock_claimed_job}
        exception = Exception("some exception")
        self.mocks['start_job_execution'].side_effect = exception
        self.runner.process_claimable_job(job)
        self.assertEqual(
            self.mocks['update_job'].call_args,
            call(job=mock_claimed_job, updates={
                'status': 'FAILED',
                'error': 'Error starting job execution: %s' % exception
            }))

class BuildJobDirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = {
            'uuid': 'abcd',
            'job_spec': {
                'pre_build_actions': MagicMock(),
            }
        }

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['process_actions'])

    def test_processes_pre_build_actions(self):
        self.runner.build_job_submission(job=self.job)
        self.assertEqual(self.runner.process_actions.call_args,
                         call(actions=self.job['job_spec']['pre_build_actions'],
                              job=self.job))

    def test_calls_job_submission_factory(self):
        self.runner.build_job_submission(job=self.job)
        self.assertEqual(
            self.runner.job_submission_factory.build_job_submission.call_args,
            call(job=self.job))

class ProcessActionsTestCase(BaseTestCase):
    def test_wraps_action_processor(self):
        actions = [MagicMock() for i in range(3)]
        job = MagicMock()
        self.runner.process_actions(actions=actions, job=job)
        self.assertEqual(
            self.runner.action_processor.process_action.call_args_list,
            [call(action=action, ctx=job) for action in actions]
        )

class StartJobExecutionTestCase(BaseTestCase):
    def test_start_job_execution(self):
        job = {'uuid': 'abcd'}
        execution_meta = self.runner.start_job_execution(job=job)
        self.assertEqual(self.runner.executing_jobs[job['uuid']], job)
        self.assertEqual(self.execution_client.start_execution.call_args,
                         call(job=job))
        self.assertEqual(execution_meta,
                         self.execution_client.start_execution.return_value)

class ProcessExecutingJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.executed_jobs = [{'uuid': 'executed_%s' % i} for i in range(3)]
        self.executing_jobs = [{'uuid': 'executing_%s' % i} for i in range(3)]
        self.runner.executing_jobs = {
            j['uuid']: j for j in self.executed_jobs + self.executing_jobs
        }
        self.mock_job_execution_states = {
            **{j['uuid']: None for j in self.executed_jobs},
            **{j['uuid']: {'executing': True} for j in self.executing_jobs}
        }
        self.mocks['get_job_execution_states'].return_value = \
                self.mock_job_execution_states
        self.runner.process_executing_jobs()

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['get_job_execution_states',
                                             'process_executed_job'])

    def test_process_executing_jobs(self):
        expected_calls = [call(job=job) for job in self.executed_jobs]
        calls = self.mocks['process_executed_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['uuid']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)

    def test_removes_executed_jobs_from_executing_jobs_registry(self):
        executed_job_uuids = [job['uuid'] for job in self.executed_jobs]
        jobs_overlap = set(self.runner.executing_jobs.keys()).intersection(
            set(executed_job_uuids))
        self.assertEqual(len(jobs_overlap), 0)

class GetJobExecutionStatesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.executing_jobs = {}
        for i in range(3):
            job_uuid = 'job_uuid_%s' % i
            self.runner.executing_jobs[job_uuid] = {'uuid': job_uuid}

    def test_get_job_execution_states(self):
        job_execution_states = self.runner.get_job_execution_states()
        expected_get_job_execution_state_calls = [
            call(job=job) for job in self.runner.executing_jobs.values()
        ]
        self.assertEqual(
            self.execution_client.get_execution_state.call_args_list,
            expected_get_job_execution_state_calls)
        expected_execution_states = {
            job['uuid']: self.execution_client.get_execution_state.return_value
            for job in self.runner.executing_jobs.values()
        }
        self.assertEqual(job_execution_states, expected_execution_states)

class ProcessExecutedJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = defaultdict(MagicMock)

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['process_actions', 'complete_job'])

    def test_adds_completed_dir_to_job_state(self):
        self.job['execution'] = {'completed_dir': 'some_dir'}
        self.runner.process_executed_job(job=self.job)
        self.assertEqual(self.job['completed_dir'],
                         self.job['execution']['completed_dir'])

    def test_calls_post_exec_actions(self):
        actions = [MagicMock() for i in range(3)]
        self.job['job_spec'] = {'post_exec_actions': actions}
        self.runner.process_executed_job(job=self.job)
        self.assertEqual(
            self.runner.process_actions.call_args,
            call(actions=self.job['job_spec']['post_exec_actions'],
                 job=self.job)
        )

    def test_calls_complete_job(self):
        self.runner.process_executed_job(job=self.job)
        self.assertEqual(self.runner.complete_job.call_args, call(job=self.job))

class CompleteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = MagicMock()
        self.runner.register_job(job=self.job)
        self.runner.complete_job(job=self.job)

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['update_job', 'unregister_job'])

    def test_updates_job(self):
        self.assertEqual(
            self.runner.update_job.call_args,
            call(job=self.job,
                 updates={
                     'status': 'COMPLETED',
                     'data': self.job.get('data')
                 }))

    def test_unregisters_job(self):
        self.assertEqual(self.runner.unregister_job.call_args,
                         call(job=self.job))

class UpdateJobTestCase(BaseTestCase):
    def test_update_job(self):
        job = {'uuid': 'abcd'}
        updates = {'pie': 'blueberry', 'meat': 'beef'}
        self.runner.update_job(job=job, updates=updates)
        self.assertEqual(self.job_client.update_jobs.call_args, 
                         call(updates_by_uuid={job['uuid']: updates}))

if __name__ == '__main__':
    unittest.main()
