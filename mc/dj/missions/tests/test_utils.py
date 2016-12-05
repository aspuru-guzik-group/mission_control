from unittest.mock import call, patch, MagicMock, DEFAULT

from django.test import TestCase

from jobs.models import Job, JobStatuses
from ..models import (Mission, Workflow, WorkflowJob, WorkflowSpec,
                      WorkflowStatuses)
from .. import utils


class WorkflowBaseTestCase(TestCase):
    def start_patchers(self, patchers):
        mocks = {key: patcher.start() for key, patcher in patchers.items()}
        return mocks

    def tearDown(self):
        if hasattr(self, 'patchers'):
            for patcher in self.patchers.values(): patcher.stop()

class StartWorkflowTestCase(WorkflowBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(
                utils, run_workflow_tick=DEFAULT,
                generate_workflow_from_spec=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.spec = WorkflowSpec.objects.create()
        self.mission = Mission.objects.create()

    def test_saves_workflow_model(self):
        workflow = utils.start_workflow(spec_key=self.spec.key,
                                        mission=self.mission)
        self.assertEqual(workflow.spec, self.spec)
        self.assertEqual(workflow.mission, self.mission)
        self.assertEqual(workflow.status, WorkflowStatuses.Running.name)
        self.assertEqual(workflow.save.call_count, 1)

    def test_calls_workflow_tick(self):
        workflow = utils.start_workflow(spec_key=self.spec.key,
                                        mission=self.mission)
        self.assertEqual(
            self.mocks['utils']['run_workflow_tick'].call_args,
            call(workflow=workflow))

class GenerateWorkflowFromSpecTestCase(WorkflowBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_engine=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.spec = WorkflowSpec.objects.create()

    def test_calls_engine_generate_workflow(self):
        workflow = utils.generate_workflow_from_spec(spec=self.spec)
        generate_workflow_fn = self.mocks['utils']['get_workflow_engine']\
                .return_value.generate_workflow_from_spec
        self.assertEqual(workflow, generate_workflow_fn.return_value)
        self.assertEqual(generate_workflow_fn.call_args, call(spec=self.spec))

class RunWorkflowTickTestCase(WorkflowBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_engine=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_calls_engine_run_workflow_tick(self):
        workflow = MagicMock()
        utils.run_workflow_tick(workflow=workflow)
        run_workflow_tick_fn = self.mocks['utils']['get_workflow_engine']\
                .return_value.run_workflow_tick
        self.assertEqual(run_workflow_tick_fn.call_args,
                         call(workflow=workflow))

class PollWorkflowJobsTestCase(WorkflowBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, run_workflow_tick=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
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
        call_args_list = self.mocks['utils']['run_workflow_tick'].call_args_list
        expected_call_args_list = [
            call(workflow=workflow_job.workflow)
            for workflow_job in self.workflow_jobs['should_finish']
        ]
        self.assertEqual(call_args_list, expected_call_args_list)
