from django.test import TestCase

from ..models import (Mission, Task, Workflow, WorkflowRunner)

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

    def test_last_finished_job(self):
        self.fail()

class WorkflowRunnerTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'name': 'test_name',
        }
        runner = WorkflowRunner.objects.create(**kwargs)
        self.assertEqual(runner.name, kwargs['name'])
        self.assertTrue(runner.uuid is not None)
        self.assertTrue(runner.created is not None)
        self.assertTrue(runner.modified is not None)
