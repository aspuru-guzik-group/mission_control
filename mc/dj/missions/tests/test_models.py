from django.test import TestCase

from jobs.models import Job
from ..models import (Mission, WorkflowStatuses, Workflow, WorkflowJob)

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

class WorkflowTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'mission': Mission.objects.create(name='mission'),
            'serialization': 'some serialization',
        }
        workflow = Workflow.objects.create(**kwargs)
        for kwarg, value in kwargs.items():
            self.assertEqual(getattr(workflow, kwarg), value)
        self.assertTrue(workflow.uuid is not None)
        self.assertTrue(workflow.created is not None)
        self.assertTrue(workflow.modified is not None)
        self.assertEqual(workflow.status, WorkflowStatuses.PENDING.name)
        self.assertTrue(hasattr(workflow, 'workflow_jobs'))
        self.assertTrue(hasattr(workflow, 'claimed'))

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
        self.assertEqual(workflow_job.meta, {})
