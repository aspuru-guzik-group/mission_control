from django.test import TestCase
from unittest.mock import call, MagicMock, patch, DEFAULT

from jobs.models import Job, JobStatuses
from ..models import Workflow, WorkflowJob
from .. import utils

class StartWorkflowTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_runner=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()
        super().tearDown()

    def test_creates_workflow(self):
        runner_key = 'some runner'
        workflow = utils.start_workflow(runner_key=runner_key)
        self.assertEqual(list(Workflow.objects.all()), [workflow])
        self.assertEqual(workflow.runner_key, runner_key)

    def test_calls_runner_tick(self):
        workflow = utils.start_workflow(runner_key='junk')
        mock_runner = self.mocks['utils']['get_workflow_runner'].return_value
        self.assertEqual(mock_runner.run_workflow_tick.call_args,
                         call(workflow=workflow))

class GetWorkflowRunnerTestCase(TestCase):
    def test_calls_loader(self):
        runner_loader = MagicMock()
        runner_key = 'some runner'
        utils.get_workflow_runner(runner_key=runner_key,
                                  runner_loader=runner_loader)
        self.assertEqual(runner_loader.load_runner.call_args, 
                         call(key=runner_key))

class PollWorkflowJobsTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(
                utils, run_workflow_job_workflow_ticks=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.workflow_jobs = self.generate_workflow_jobs()

    def generate_workflow_jobs(self):
        keyed_jobs = {
            'should_finish': self.generate_finished_jobs(),
            'should_not_finish': self.generate_unfinished_jobs(),
        }
        keyed_workflow_jobs = {
            key: [
                WorkflowJob.objects.create(
                    workflow=Workflow.objects.create(),
                    job=job
                ) for job in jobs_for_key
            ]
            for key, jobs_for_key in keyed_jobs.items()
        }
        keyed_workflow_jobs['already_finished'] = [
            WorkflowJob.objects.create(
                workflow=Workflow.objects.create(),
                job=Job.objects.create(),
                finished=True
            )
        ]
        return keyed_workflow_jobs

    def generate_finished_jobs(self):
        finished_statuses = [JobStatuses.Completed.name]
        return [Job.objects.create(status=status)
                for status in finished_statuses]

    def generate_unfinished_jobs(self):
        unfinished_statuses = [None, JobStatuses.Claimed.name,
                               JobStatuses.Pending.name]
        return [Job.objects.create(status=status)
                for status in unfinished_statuses]

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()
        super().tearDown()

    def test_marks_jobs_as_finished(self):
        utils.poll_workflow_jobs()
        finished_workflow_jobs = WorkflowJob.objects.filter(finished=True)
        expected_finished_workflow_jobs = [
            *self.workflow_jobs['should_finish'],
            *self.workflow_jobs['already_finished']
        ]
        self.assertEqual(set(finished_workflow_jobs),
                         set(expected_finished_workflow_jobs))
        unfinished_workflow_jobs = WorkflowJob.objects.filter(finished=False)
        expected_unfinished_workflow_jobs = self.workflow_jobs.get(
            'should_not_finish')
        self.assertEqual(set(unfinished_workflow_jobs),
                         set(expected_unfinished_workflow_jobs))

    def test_runs_workflow_ticks_for_jobs_that_finished(self):
        utils.poll_workflow_jobs()
        mock = self.mocks['utils']['run_workflow_job_workflow_ticks']
        workflow_jobs_kwarg = mock.call_args[1]['workflow_jobs']
        expected_workflow_jobs_kwarg = self.workflow_jobs['should_finish']
        self.assertEqual(set(workflow_jobs_kwarg),
                         set(expected_workflow_jobs_kwarg))

class RunWorkflowJobsWorkflowTicksTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_runner=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.workflow_jobs = self.generate_workflow_jobs()

    def generate_workflow_jobs(self):
        workflow_jobs = [
            WorkflowJob(
                workflow=Workflow.objects.create(runner_key=(i % 2)),
                job=Job.objects.create()
            )
            for i in range(6)
        ]
        return workflow_jobs

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()
        super().tearDown()

    def test_runs_workflow_ticks(self):
        utils.run_workflow_job_workflow_ticks(workflow_jobs=self.workflow_jobs)
        mock_runner = self.mocks['utils']['get_workflow_runner'].return_value
        call_args_list = mock_runner.run_workflow_tick.call_args_list
        expected_call_args_list = [call(workflow=workflow_job.workflow)
                                   for workflow_job in self.workflow_jobs]
        self.assertEqual(
            sorted(call_args_list, key=lambda c: c[1]['workflow'].uuid),
            sorted(expected_call_args_list, key=lambda c: c[2]['workflow'].uuid)
        )

class SyncWorkflowRunnerRegistryTestCase(TestCase):
    def test_registers_new_runners(self):
        utils.sync_workflow_runner_registry()
        runners = WorkflowRunner.objects.all()
        expected_runners = None
        self.assertEqual(set(runners), set(expected_runners))

    def test_removes_deleted_runners(self):
        self.fail()

    def test_handles_registration_errors(self):
        self.fail()
