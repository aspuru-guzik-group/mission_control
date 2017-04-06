from collections import defaultdict
import logging
import unittest
from unittest.mock import call, MagicMock, patch
from ..base_job_runner import BaseJobRunner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = defaultdict(MagicMock)
        self.task_handler = MagicMock()
        self.job_client = MagicMock()
        self.job_submission_factory = MagicMock()
        self.runner = self.generate_base_job_runner()
        self.runner_methods_to_patch = []

    def generate_base_job_runner(self, **runner_kwargs):
        runner_kwargs = runner_kwargs or {}
        default_kwargs = {
            'task_handler': self.task_handler,
            'job_client': self.job_client,
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
    def setUp(self):
        super().setUp()
        for method_name in ['fetch_claimable_jobs', 'process_claimable_job',
                            'process_running_jobs']:
            setattr(self.runner, method_name, MagicMock())

    def test_processes_running_jobs(self):
        self.runner.tick()
        self.assertEqual(self.runner.process_running_jobs.call_count, 1)

    def test_fetches_claimable_job(self):
        self.runner.tick()
        self.assertEqual(self.runner.fetch_claimable_jobs.call_count, 1)

    def test_processes_claimable_jobs(self):
        self.runner.max_running_jobs = 3
        extra_jobs = 2
        jobs = [
            {'uuid': i}
            for i in range(self.runner.max_running_jobs + extra_jobs)
        ]
        self.runner.fetch_claimable_jobs.return_value = jobs
        self.runner.tick()
        expected_call_args = [call(job=job) for job in jobs[:-1 * extra_jobs]]
        self.assertTrue(self.runner.process_claimable_job.call_args_list,
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
        self.runner.complete_job = MagicMock()
        self.runner.process_running_jobs()

    def test_completes_expected_jobs(self):
        expected_calls = [call(job=job) for job in self.completed_jobs]
        def _calls_by_job_uuid(calls):
            return {_call[-1]['job']['uuid']: _call for _call in calls}
        self.assertEqual(
            _calls_by_job_uuid(self.runner.complete_job.call_args_list),
            _calls_by_job_uuid(expected_calls)
        )

class TickJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.get_task_runner = MagicMock()
        self.expected_task_runner = self.runner.get_task_runner.return_value

    def test_dispatches_to_task_runner_tick_tasks(self):
        self.runner.tick_job(job=self.job)
        self.assertEqual(self.runner.get_task_runner.call_args,
                         call(job=self.job))
        self.assertEqual(self.expected_task_runner.tick_tasks.call_args, call())

    def test_fails_job_on_execption(self):
        exception = Exception('some exception')
        self.expected_task_runner.tick_tasks.side_effect = exception
        self.runner.fail_job = MagicMock()
        self.runner.tick_job(job=self.job)
        self.assertEqual(
            self.runner.fail_job.call_args,
            call(job=self.job, error=self.runner.stringify_exception(exception))
        )

class CompleteJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = defaultdict(MagicMock)
        self.job['execution'] = defaultdict(MagicMock)
        self.runner.register_job(job=self.job)
        self.runner.update_job = MagicMock()
        self.runner.unregister_job = MagicMock()

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
        self.job['status'] = 'FAILED'
        self.job['error'] = 'some error'
        self.runner.fail_job = MagicMock()
        self.runner.complete_job(job=self.job)
        self.assertEqual(self.runner.fail_job.call_args,
                         call(job=self.job, error=self.job['error']))

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
    def setUp(self):
        super().setUp()
        for method_name in ['start_job', 'update_job']:
            setattr(self.runner, method_name, MagicMock())

    def tearDown(self):
        logging.disable(logging.NOTSET)
        super().tearDown()

    def test_claimable_job(self):
        job = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {job['uuid']: job}
        self.runner.process_claimable_job(job=job)
        self.assertEqual(self.runner.start_job.call_args, call(job=job))

    def test_unclaimable_job(self):
        job = {'uuid': 'abcd'}
        self.job_client.claim_jobs.return_value = {job['uuid']: None}
        self.runner.process_claimable_job(job)
        self.assertEqual(self.runner.start_job.call_args, None)

class StartJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.get_default_job_tasks = MagicMock()
        self.runner = self.generate_base_job_runner(
            get_default_job_tasks=self.get_default_job_tasks)
        self.runner.tick_job = MagicMock()

    def test_ticks_job_tasks(self):
        self.runner.start_job(job=self.job)
        self.assertEqual(self.runner.tick_job.call_args, call(job=self.job))

    def test_fallsback_to_default_tasks(self):
        self.runner.start_job(job=self.job)
        self.assertEqual(self.get_default_job_tasks.call_args,
                         call(job=self.job))
        self.assertEqual(self.job['tasks'], self.get_default_job_tasks())

class GetTaskContextTestCase(BaseTestCase):
    def _do_get_task_context(self):
        return self.runner.get_task_context(job=self.job)

    def test_context_has_job(self):
        task_context = self._do_get_task_context()
        self.assertEqual(task_context['job'], self.job)

if __name__ == '__main__':
    unittest.main()
