from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock, patch
from .. import base_job_runner


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
        return base_job_runner.BaseJobRunner(
            **{**default_kwargs, **runner_kwargs})

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
        for method_name in ['claim_jobs', 'process_job',
                            'process_running_jobs']:
            setattr(self.runner, method_name, MagicMock())

    def test_processes_running_jobs(self):
        self.runner.tick()
        self.assertEqual(self.runner.process_running_jobs.call_count, 1)

    def test_claims_jobs(self):
        self.runner.tick()
        self.assertEqual(self.runner.claim_jobs.call_args,
                         call(params={'limit': self.runner.max_running_jobs}))

class ProcessRunningJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.running_jobs = [{'uuid': 'running_%s' % i} for i in range(3)]
        self.runner.running_jobs = {job['uuid']: job
                                    for job in self.running_jobs}
        for attr in ['tick_job', 'update_jobs']:
            setattr(self.runner, attr, MagicMock())
        def mock_job_is_running(job=None): return ('running' in job['uuid'])
        self.runner.job_is_running = MagicMock(side_effect=mock_job_is_running)
        self.runner.process_running_jobs()

    def test_ticks_running_jobs(self):
        expected_calls = [call(job=job) for job in self.running_jobs]
        self.assertEqual(
            self._calls_by_job_uuid(self.runner.tick_job.call_args_list),
            self._calls_by_job_uuid(expected_calls)
        )

    def _calls_by_job_uuid(self, calls):
        return {_call[-1]['job']['uuid']: _call for _call in calls}

    def test_updates_running_jobs(self):
        self.assertEqual(
            set([
                j['uuid'] for j in self.runner.update_jobs.call_args[1]['jobs']
            ]),
            set([j['uuid'] for j in self.running_jobs])
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

    @patch.object(base_job_runner, 'traceback')
    def test_fails_job_on_exception(self, mock_traceback):
        exception = Exception('some exception')
        self.expected_task_runner.tick_tasks.side_effect = exception
        self.runner.fail_job = MagicMock()
        self.runner.tick_job(job=self.job)
        self.assertEqual(
            self.runner.fail_job.call_args,
            call(job=self.job, error=mock_traceback.format_exc.return_value))

class UpdateJobsTestCase(BaseTestCase):
    def test_update_job(self):
        jobs = [MagicMock() for i in range(3)]
        self.runner.update_jobs(jobs=jobs)
        self.assertEqual(
            self.job_client.update_jobs.call_args, 
            call(updates_by_uuid={
                job['uuid']: {attr: job.get(attr)
                              for attr in ['data', 'status', 'error']}
                for job in jobs
            })
        )

class ClaimJobsTestCase(BaseTestCase):
    def test_fetch_claimable_jobs(self):
        params = MagicMock()
        self.runner.claim_jobs(params=params)
        self.assertEqual(self.job_client.claim_jobs.call_args,
                         call(params=params))

class StartJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.get_default_job_tasks = MagicMock()
        self.runner = self.generate_base_job_runner(
            get_default_job_tasks=self.get_default_job_tasks)
        self.runner.tick_job = MagicMock()

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
