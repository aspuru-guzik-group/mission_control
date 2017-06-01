from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import post_chemthing_actions_task_handler

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task_handler = post_chemthing_actions_task_handler.TaskHandler()
        self.task = defaultdict(MagicMock)
        self.task_context = MagicMock()

@patch.object(post_chemthing_actions_task_handler, '_a2g2_db_utils')
class InitialTickTestCase(BaseTestCase):
    def test_posts_chemthing_actions_via_client(self, mock_a2g2_db_utils):
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)
        self.assertEqual(mock_a2g2_db_utils.generate_a2g2_client.call_args,
                         call(cfg=self.task_context['flow_ctx']['cfg']))
        expected_client = mock_a2g2_db_utils.generate_a2g2_client.return_value
        expected_chemthing_actions = \
                self.task['task_params']['chemthing_actions']
        self.assertEqual(expected_client.post_chemthing_actions.call_args,
                         call(chemthing_actions=expected_chemthing_actions))

    def test_sets_task_status_to_completed(self, mock_a2g2_db_utils):
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)
        self.assertEqual(self.task['status'], 'COMPLETED')
