from collections import defaultdict
import logging
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from ..base_job_runner import BaseJobRunner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = defaultdict(MagicMock)
        self.task_runner = MagicMock()
        self.execution_client = MagicMock()
        self.job_client = MagicMock()
        self.job_submission_factory = MagicMock()
        self.runner = self.generate_base_job_runner()
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

    def generate_base_job_runner(self, runner_kwargs=None):
        runner_kwargs = runner_kwargs or {}
        default_kwargs = {
            'task_runner': self.task_runner,
            'execution_client': self.execution_client,
            'job_client': self.job_client,
            'job_submission_factory': self.job_submission_factory,
        }
        return BaseJobRunner(**{**default_kwargs, **runner_kwargs})

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
                                             'process_running_jobs'])

    def test_processes_running_jobs(self):
        self.runner.tick()
        self.assertEqual(self.mocks['process_running_jobs'].call_count, 1)

    def test_fetches_claimable_job(self):
        self.runner.tick()
        self.assertTrue(self.mocks['fetch_claimable_jobs'].called)

    def test_processes_claimable_jobs(self):
        self.runner.max_running_jobs = 3
        extra_jobs = 2
        jobs = [
            {'uuid': i}
            for i in range(self.runner.max_running_jobs + extra_jobs)
        ]
        self.mocks['fetch_claimable_jobs'].return_value = jobs
        self.runner.tick()
        expected_call_args = [call(job=job)
                              for job in jobs[:-1 * extra_jobs]]
        self.assertTrue(self.mocks['process_claimable_job'].call_args_list,
                        expected_call_args)

class ProcessRunningJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.completed_jobs = [{'uuid': 'completed_%s' % i} for i in range(3)]
        self.running_jobs = [{'uuid': 'running_%s' % i} for i in range(3)]
        self.runner.running_jobs = {
            j['uuid']: j for j in self.completed_jobs + self.running_jobs
        }
        def mock_job_is_running(job=None): return ('running' in job['uuid'])
        self.runner.job_is_running = mock_job_is_running
        self.runner.process_running_jobs()

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['complete_job'])

    def test_completes_expected_jobs(self):
        expected_calls = [call(job=job) for job in self.completed_jobs]
        calls = self.mocks['complete_job'].mock_calls
        calls_sort_key_fn = lambda c: c[2]['job']['uuid']
        sorted_calls = sorted(calls, key=calls_sort_key_fn)
        sorted_expected_calls = sorted(expected_calls, key=calls_sort_key_fn)
        self.assertEqual(sorted_calls, sorted_expected_calls)

class CompleteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = defaultdict(MagicMock)
        self.job['execution'] = defaultdict(MagicMock)
        self.runner.register_job(job=self.job)

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['update_job', 'unregister_job'])

    def test_updates_successful_job_and_unregisters(self):
        self.job['execution']['result'] = {'result': 'COMPLETED'}
        self.job['data'] = 'some data'
        self.runner.complete_job(job=self.job)
        self.assertEqual(
            self.runner.update_job.call_args,
            call(job=self.job,
                 updates={
                     'status': 'COMPLETED',
                     'data': self.job['data']
                 }))
        self.assertEqual(self.runner.unregister_job.call_args,
                         call(job=self.job))

    def test_fails_failed_job(self):
        error = 'some error'
        self.job['execution']['result'] = {'result': 'FAILED', 'error': error}
        self.runner.fail_job = MagicMock()
        self.runner.complete_job(job=self.job)
        self.assertEqual(self.runner.fail_job.call_args,
                         call(job=self.job, error=error))

class UpdateJobTestCase(BaseTestCase):
    def test_update_job(self):
        job = {'uuid': 'abcd'}
        updates = {'pie': 'blueberry', 'meat': 'beef'}
        self.runner.update_job(job=job, updates=updates)
        self.assertEqual(self.job_client.update_jobs.call_args, 
                         call(updates_by_uuid={job['uuid']: updates}))

class FetchClaimableJobsTestCase(BaseTestCase):
    def test_fetch_claimable_jobs(self):
        self.runner.fetch_claimable_jobs()
        self.assertEqual(
            self.job_client.fetch_claimable_jobs.call_count, 1)

class ProcessClaimableJobTestCase(BaseTestCase):
    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['build_job_submission',
                                             'start_job',
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
        self.assertEqual(self.mocks['build_job_submission'].call_count, 1)
        self.assertEqual(job['submission'], expected_submission_meta)

    def test_unclaimable_job(self):
        job = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {
            job['uuid']: None}
        self.runner.process_claimable_job(job)
        self.assertEqual(self.mocks['start_job'].call_count, 0)

class BuildJobSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = {'uuid': 'abcd'}

    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['generate_job_dir'])

    def test_generates_job_dir(self):
        self.runner.build_job_submission(job=self.job)
        self.assertEqual(self.job['job_dir'],
                         self.runner.generate_job_dir.return_value)

    def test_calls_job_submission_factory(self):
        self.runner.build_job_submission(job=self.job)
        self.assertEqual(
            self.runner.job_submission_factory.build_job_submission.call_args,
            call(job=self.job, output_dir=self.job['job_dir'])
        )

class StartJobTestCase(BaseTestCase):
    def decorate_runner_methods_to_patch(self):
        self.runner_methods_to_patch.extend(['tick_job_tasks'])

    def test_ticks_job_tasks(self):
        self.runner.start_job(job=self.job)
        self.assertEqual(self.mocks['tick_job_tasks'].call_args,
                         call(job=self.job))

    def test_fallsback_to_default_tasks(self):
        self.runner.default_tasks = MagicMock()
        self.runner.start_job(job=self.job)
        self.assertEqual(self.job['tasks'], self.runner.default_tasks)

class TickJobTasksTestCase(BaseTestCase):
    def test_continues_until_first_noncompleting_task(self):
        completing_tasks = [{'key': i} for i in range(3)]
        noncompleting_tasks = [{'key': i} for i in range(3)]
        tasks = completing_tasks + noncompleting_tasks
        self.job['tasks'] = tasks
        def mock_tick_task(task=None, job=None):
            if task in completing_tasks: task['status'] = 'COMPLETED'
            else: task['status'] = 'RUNNING'
        self.runner.tick_task = MagicMock(side_effect=mock_tick_task)
        self.runner.tick_job_tasks(job=self.job)
        expected_tick_task_call_args_list = (
            [call(task=task, job=self.job) for task in completing_tasks] +
            [call(task=noncompleting_tasks[0], job=self.job)]
        )
        self.assertEqual(self.runner.tick_task.call_args_list,
                         expected_tick_task_call_args_list)

class TickTaskTestCase(BaseTestCase):
    def test_dispatches_to_task_runner(self):
        task = MagicMock()
        self.runner.tick_task(task=task, job=self.job)
        self.assertEqual(self.task_runner.tick_task.call_args,
                         call(task=task, job=self.job))

if __name__ == '__main__':
    unittest.main()
