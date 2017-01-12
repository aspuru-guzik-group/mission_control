from django.test import TestCase

from jobs.models import Job
from ..models import (Mission, FlowStatuses, Flow, FlowJob)

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

class FlowTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'mission': Mission.objects.create(name='mission'),
            'serialization': 'some serialization',
            'spec': 'some spec',
        }
        flow = Flow.objects.create(**kwargs)
        for kwarg, value in kwargs.items():
            self.assertEqual(getattr(flow, kwarg), value)
        self.assertTrue(flow.uuid is not None)
        self.assertTrue(flow.created is not None)
        self.assertTrue(flow.modified is not None)
        self.assertEqual(flow.status, FlowStatuses.PENDING.name)
        self.assertTrue(hasattr(flow, 'flow_jobs'))
        self.assertTrue(hasattr(flow, 'claimed'))

class FlowJobTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'flow': Flow.objects.create(),
            'job': Job.objects.create(),
        }
        flow_job = FlowJob.objects.create(**kwargs)
        self.assertTrue(flow_job.uuid is not None)
        self.assertTrue(flow_job.created is not None)
        self.assertTrue(flow_job.modified is not None)
        self.assertTrue(flow_job.flow is kwargs['flow'])
        self.assertTrue(flow_job.job is kwargs['job'])
        self.assertEqual(flow_job.finished, False)
        self.assertEqual(flow_job.meta, {})
