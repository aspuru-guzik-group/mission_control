from collections import defaultdict
import unittest
from unittest.mock import MagicMock

from .. import base_flow_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = defaultdict(MagicMock, **{'data': {}})
        self.flow_engine = MagicMock()
        self.task_ctx = self.generate_task_ctx(flow_engine=self.flow_engine)
        self.task_handler = base_flow_task_handler.BaseFlowTaskHandler(
            task_ctx=self.task_ctx)

    def generate_task_ctx(self, **kwargs):
        task_ctx = defaultdict(MagicMock, **{'task': self.task, **kwargs})
        return task_ctx

class HandleFlowStatusTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = MagicMock()

    def _handle(self):
        self.task_handler.handle_flow_status(flow=self.flow)

    def test_sets_status_and_data_for_completed_flow(self):
        self.flow.status = 'COMPLETED'
        self._handle()
        self.assertEqual(self.task_handler.task['status'], 'COMPLETED')
        self.assertEqual(self.task_handler.task['data']['flow_data'],
                         self.flow.data)

    def test_sets_errors_for_failed_flow_and_raises(self):
        self.flow.status = 'FAILED'
        expected_msg = str(self.flow.data.get('errors', '<unknown>'))
        with self.assertRaises(Exception, msg=expected_msg):
            self._handle()

if __name__ == '__main__':
    unittest.main()
