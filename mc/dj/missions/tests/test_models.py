from django.test import TestCase

from jobs.models import Job
from ..models import Mission, Task, Workflow, WorkflowJob, WorkflowRunner

class MissionTestCase(TestCase):
    def test_has_expected_fields(self):
        mission_kwargs = {
            'name': 'test_name',
        }
        mission = Mission.objects.create(**mission_kwargs)
        self.assertEqual(mission.name, mission_kwargs['name'])
        self.assertTrue(mission.uuid is not None)
        self.assertTrue(mission.created is not None)
        self.assertTrue(mission.modified is not None)

class TaskTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'name': 'name',
            'mission': Mission.objects.create(name='name'),
        }
        task = Task.objects.create(**kwargs)
        self.assertEqual(task.mission, kwargs['mission'])
        self.assertTrue(task.uuid is not None)
        self.assertTrue(task.created is not None)
        self.assertTrue(task.modified is not None)

class WorkflowTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'runner_key': 'runner_key',
            'mission': Mission.objects.create(name='mission'),
        }
        workflow = Workflow.objects.create(**kwargs)
        self.assertEqual(workflow.mission, kwargs['mission'])
        self.assertTrue(workflow.uuid is not None)
        self.assertTrue(workflow.created is not None)
        self.assertTrue(workflow.modified is not None)
        self.assertTrue(workflow.jobs is not None)

    def test_last_finished_job(self):
        workflow = Workflow.objects.create()
        for i in range(3):
            WorkflowJob.objects.create(
                workflow=workflow,
                job=Job.objects.create(),
                finished=True,
                modified=i)
        last_finished_job = workflow.last_finished_job
        expected_last_finished_job = workflow.workflow_jobs.order_by(
            '-modified').first().job
        self.assertEqual(last_finished_job, expected_last_finished_job)

class WorkflowJobTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'workflow': Workflow.objects.create(),
            'job': Job.objects.create(),
        }
        workflow_job = WorkflowJob.objects.create(**kwargs)
        self.assertTrue(workflow_job.uuid is not None)
        self.assertTrue(workflow_job.created is not None)
        self.assertTrue(workflow_job.modified is not None)
        self.assertTrue(workflow_job.workflow is kwargs['workflow'])
        self.assertTrue(workflow_job.job is kwargs['job'])
        self.assertEqual(workflow_job.finished, False)

class WorkflowRunnerTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'key': 'some key',
            'path': 'some path',
            'label': 'some label',
        }
        workflow_runner = WorkflowRunner.objects.create(**kwargs)
        self.assertTrue(workflow_runner.uuid is not None)
        self.assertTrue(workflow_runner.created is not None)
        self.assertTrue(workflow_runner.modified is not None)
        for attr, value in kwargs.items():
            self.assertEqual(getattr(workflow_runner, attr), value)
