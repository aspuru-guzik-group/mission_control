import collections
from unittest.mock import call, patch, MagicMock, DEFAULT

from django.test import TestCase

from jobs.models import Job, JobStatuses
from ..models import (Mission, Workflow, WorkflowJob, WorkflowSpec,
                      WorkflowStatuses)
from .. import utils


class BaseTestCase(TestCase):
    def start_patchers(self, patchers):
        mocks = {key: patcher.start() for key, patcher in patchers.items()}
        return mocks

    def tearDown(self):
        if hasattr(self, 'patchers'):
            for patcher in self.patchers.values(): patcher.stop()

class StartWorkflowTestCase(BaseTestCase):
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

class GenerateWorkflowFromSpecTestCase(BaseTestCase):
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

class RunWorkflowTickTestCase(BaseTestCase):
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

class PollWorkflowJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils,
                                    process_finished_workflow_jobs=DEFAULT)
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

    def test_processes_finished_jobs(self):
        utils.poll_workflow_jobs()
        self.assertEqual(
            self.mocks['utils']['process_finished_workflow_jobs'].call_args,
            call(workflow_jobs=self.workflow_jobs['should_finish']))

class ProcessFinishedWorkflowJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_engine=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.engine = self.mocks['utils']['get_workflow_engine'].return_value
        self.workflow_jobs = self.generate_workflow_jobs()
        self.grouped_workflow_jobs = collections.defaultdict(list)
        for workflow_job in self.workflow_jobs:
            self.grouped_workflow_jobs[workflow_job.workflow].append(
                workflow_job)

    def generate_workflow_jobs(self):
        workflows = [Workflow.objects.create() for i in range(3)]
        workflow_jobs = []
        for workflow in workflows:
            workflow_jobs.extend([WorkflowJob.objects.create(
                workflow=workflow,
                job=Job.objects.create(),
            ) for i in range(3)])
        return workflow_jobs

    def test_calls_engine_process_finished_workflow_jobs(self):
        utils.process_finished_workflow_jobs(workflow_jobs=self.workflow_jobs)
        call_args_list = self.engine.process_finished_workflow_jobs\
                .call_args_list
        expected_call_args_list = [
            call(workflow=workflow, workflow_jobs=workflow_jobs)
            for workflow, workflow_jobs in self.grouped_workflow_jobs.items()
        ]
        self.assertEqual(call_args_list, expected_call_args_list)

    def test_runs_workflow_tick(self):
        utils.process_finished_workflow_jobs(workflow_jobs=self.workflow_jobs)
        call_args_list = self.engine.run_workflow_tick.call_args_list
        expected_call_args_list = [
            call(workflow=workflow)
            for workflow in self.grouped_workflow_jobs.keys()
        ]
        self.assertEqual(call_args_list, expected_call_args_list)

class CreateSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_engine=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mock_engine = self.mocks['utils']['get_workflow_engine']\
                .return_value
        self.mock_engine.normalize_spec_serialization.return_value = (
            'normalized')
        self.serialization = 'some serialization'
        self.key = 'some key'

    def test_returns_expected_spec(self):
        expected_serialization = self.mock_engine.normalize_spec_serialization\
                .return_value
        created_spec = utils.create_spec(key=self.key,
                                         serialization=expected_serialization)
        self.assertEqual(created_spec.key, self.key)

    def test_persists_created_spec(self):
        created_spec = utils.create_spec(serialization=self.serialization)
        specs_from_db = WorkflowSpec.objects.all()
        self.assertEqual(set(specs_from_db), set([created_spec]))

