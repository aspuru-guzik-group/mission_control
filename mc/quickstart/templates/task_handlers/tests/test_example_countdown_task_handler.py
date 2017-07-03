import unittest

from ..example_countdown_task_handler import TaskHandler


class ExampleCountdownTaskHandlerTestCase(unittest.TestCase):
    def test_my_countdown_task_handler(self):
        countdown_start = 3
        task = {'key': 'test_key', 'data': {},
                'task_params': {'countdown_start': countdown_start}}
        expected_num_ticks = countdown_start + 1
        tick_counter = 0
        while task.get('status') not in {'COMPLETED', 'FAILED'}:
            TaskHandler.tick_task(task_ctx={'task': task})
            tick_counter += 1
            if tick_counter > expected_num_ticks:
                raise Exception("Exceeded expected number of ticks")
        self.assertEqual(task['data']['result'], 'Lift off')
        self.assertEqual(task['status'], 'COMPLETED')
