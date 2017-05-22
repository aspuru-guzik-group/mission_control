from django.test import TestCase

from ..constants import Statuses
from ..models import (Mission, Flow, Job, Queue)


class MissionTestCase(TestCase):
    def test_has_expected_fields(self):
        mission_kwargs = {
            'label': 'test_label',
        }
        mission = Mission.objects.create(**mission_kwargs)
        self.assertEqual(mission.label, mission_kwargs['label'])
        self.assertTrue(mission.key is not None)
        self.assertTrue(mission.created is not None)
        self.assertTrue(mission.modified is not None)

class FlowTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {'serialization': 'some serialization'}
        flow = Flow.objects.create(**kwargs)
        for kwarg, value in kwargs.items():
            self.assertEqual(getattr(flow, kwarg), value)
        self.assertTrue(flow.key is not None)
        self.assertTrue(flow.created is not None)
        self.assertTrue(flow.modified is not None)
        self.assertEqual(flow.status, Statuses.PENDING)
        self.assertTrue(hasattr(flow, 'claimed'))

class JobTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {'label': 'test_label'}
        job = Job.objects.create(**kwargs)
        self.assertEqual(job.label, kwargs['label'])
        self.assertEqual(job.status, Statuses.PENDING)
        expected_attrs = ['key', 'created', 'modified', 'claimed',
                          'serialization']
        for attr in expected_attrs: self.assertTrue(hasattr(job, attr))

class QueueTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'label': 'test_label',
            'queue_spec': {'key1': 'value1'}
        }
        queue = Queue.objects.create(**kwargs)
        self.assertEqual(queue.label, kwargs['label'])
        self.assertEqual(queue.queue_spec, kwargs['queue_spec'])
        expected_attrs = ['key', 'created', 'modified']
        for attr in expected_attrs:
            self.assertTrue(hasattr(queue, attr))
