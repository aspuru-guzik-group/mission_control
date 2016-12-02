from unittest.mock import call, patch, MagicMock, DEFAULT

from django.test import TestCase

from jobs.models import Job, JobStatuses
from ..models import Workflow, WorkflowJob, WorkflowSpec, WorkflowStatuses
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
                utils, run_workflow_tick=DEFAULT, serialize_workflow=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mocks['utils']['serialize_workflow'].return_value = 'serialization'
        self.spec = WorkflowSpec.objects.create(key='some key')

    def test_creates_workflow_model(self):
        workflow = utils.start_workflow(spec_key=self.spec.key)
        self.assertEqual(list(Workflow.objects.all()), [workflow])
        self.assertEqual(workflow.spec, self.spec)
        self.assertEqual(workflow.serialization,
                         self.mocks['utils']['serialize_workflow'].return_value)
        self.assertEqual(workflow.status, WorkflowStatuses.Running.name)

    def test_calls_workflow_tick(self):
        workflow = utils.start_workflow(spec_key=self.spec.key)
        self.assertEqual(
            self.mocks['utils']['run_workflow_tick'].call_args,
            call(workflow=workflow))

class SerializeWorkflowTestCase(WorkflowBaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_engine=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_calls_engine_serialize_workflow(self):
        workflow = MagicMock()
        serialization = utils.serialize_workflow(workflow=workflow)
        serialize_workflow_fn = self.mocks['utils']['get_workflow_engine']\
                .return_value.serialize_workflow
        self.assertEqual(serialization, serialize_workflow_fn.return_value)
        self.assertEqual(serialize_workflow_fn.call_args,
                         call(workflow=workflow))

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
