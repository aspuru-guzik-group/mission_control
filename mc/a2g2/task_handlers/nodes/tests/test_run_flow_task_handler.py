from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock
from uuid import uuid4

from ..run_flow_task_handler import RunFlowTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = MagicMock()
        self.flow_ctx = {
            'create_flow': MagicMock(),
            'get_flow': MagicMock()
        }
        self.task_handler = RunFlowTaskHandler()

    def generate_flow(self, uuid=None, status='PENDING', **flow_state):
        flow = defaultdict(MagicMock)
        if not uuid: uuid = str(uuid4())
        flow.update({'uuid': uuid, 'status': status, **flow_state})
        return flow

    def generate_task_context(self, **kwargs):
        task_context = {
            'task': self.task,
            'flow_ctx': self.flow_ctx,
            **kwargs
        }
        return task_context

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.initial_task = {
            'data': {},
            'task_params': {
                'flow_spec': 'some flow_spec'
            }
        }
        self.task = {**self.initial_task}
        self.task_handler.initial_tick(
            task=self.task,
            task_context=self.generate_task_context()
        )

    def test_initial_tick_creates_flow(self):
        self.assertEqual(
            self.flow_ctx['create_flow'].call_args,
            call(flow_kwargs={
                'flow_spec': self.initial_task['task_params']['flow_spec'],
            })
        )

    def test_has_flow_uuid(self):
        self.assertEqual(self.task['data']['flow_uuid'],
                         self.flow_ctx['create_flow'].return_value['uuid'])

class IntermediateTickMixin(object):
    def do_intermediate_tick(self, flow_state=None):
        if not flow_state: flow_state = {}
        self.flow = self.generate_flow(**flow_state)
        self.flow_ctx['get_flow'].return_value = self.flow
        self.initial_task = {
            'data': {'flow_uuid': self.flow['uuid']},
            'task_params': {'flow_spec': 'some flow spec'},
            'status': 'some_status'
        }
        self.task = {**self.initial_task}
        self.task_handler.intermediate_tick(
            task=self.task,
            task_context=self.generate_task_context()
        )

class IncompleteFlowTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.do_intermediate_tick(flow_state={'status': 'PENDING'})

    def test_has_expected_state(self):
        self.assertEqual(self.task['status'], self.initial_task['status'])
        self.assertEqual(self.task['data'], self.initial_task['data'])

class CompletedFlowTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.do_intermediate_tick(flow_state={'status': 'COMPLETED'})

    def test_has_expected_status(self):
        self.assertEqual(self.task['status'], 'COMPLETED')

class FailedFlowTestCase(BaseTestCase, IntermediateTickMixin):
    def test_throws_exception(self):
        with self.assertRaises(Exception) as ctx:
            error = 'some error'
            self.do_intermediate_tick(flow_state={
                'status': 'FAILED',
                'data': {'error': error},
            })
            self.assertTrue(error in ctx.exception)

if __name__ == '__main__':
    unittest.main()
