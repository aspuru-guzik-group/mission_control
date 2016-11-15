from django.test import TestCase

from .models import (Mission, Task)

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
        task_kwargs = {
            'name': 'name',
            'mission': Mission.objects.create(name='name'),
        }
        task = Task.objects.create(**task_kwargs)
        self.assertEqual(task.mission, task_kwargs['mission'])
        self.assertTrue(task.uuid is not None)
        self.assertTrue(task.created is not None)
        self.assertTrue(task.modified is not None)
