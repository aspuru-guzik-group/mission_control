from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from .. import collect_chemthing_actions_task_handler

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task_handler = collect_chemthing_actions_task_handler.TaskHandler()
        self.task = defaultdict(MagicMock, **{
            'data': defaultdict(MagicMock)
        })
        self.task_context = MagicMock()

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.collect_chemthing_actions = MagicMock()
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_collects_chemthing_actions(self):
        self.assertEqual(self.task_handler.collect_chemthing_actions.call_args,
                         call(task_context=self.task_context))
        self.assertEqual(
            self.task['data']['chemthing_actions'],
            self.task_handler.collect_chemthing_actions.return_value)

    def test_sets_task_status_to_completed(self):
        self.assertEqual(self.task['status'], 'COMPLETED')

class CollectChemThingActionsTestCase(BaseTestCase):
    def test_scans_node_tasks(self):
        chemthing_actions = self.task_handler.collect_chemthing_actions(
            task_context=self.task_context)
        expected_chemthing_actions = []
        for task in self.task_context['node']['node_tasks']:
            expected_chemthing_actions.extend(task['data'].get(
                'chemthing_actions', []))
        self.assertEqual(chemthing_actions, expected_chemthing_actions)
