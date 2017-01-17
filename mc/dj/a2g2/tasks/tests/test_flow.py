import json
import unittest
from unittest.mock import call, MagicMock
from uuid import uuid4
from django.test import TestCase

from ..flow import FlowTask


class BaseTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.ctx = {
            'create_flow': MagicMock(),
            'get_flow': MagicMock()
        }

    def generate_task(self, **task_kwargs):
        return FlowTask(**task_kwargs)

    def generate_flow(self, uuid=None, status='PENDING', **flow_state):
        if not uuid: uuid = str(uuid4())
        return {'uuid': uuid, 'status': status, **flow_state}

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.initial_state = {
            'input': {'flow_spec': 'some flow_spec'}
        }
        self.task = self.generate_task(**self.initial_state)
        self.task.tick(ctx=self.ctx)

    def test_initial_tick_creates_flow(self):
        expected_call_args = call(flow={
            'spec': json.dumps(self.initial_state['input']['flow_spec']),
        })
        self.assertEqual(self.ctx['create_flow'].call_args, expected_call_args)

    def test_has_expected_data(self):
        expected_data = {
            'flow_uuid': self.ctx['create_flow'].return_value['uuid'],
            'ticks': 1
        }
        self.assertEqual(self.task.data, expected_data)

    def test_has_expected_status(self):
        self.assertEqual(self.task.status, 'RUNNING')

class IntermediateTickMixin(object):
    def setup_for_intermediate_tick(self, flow_state=None):
        if not flow_state: flow_state = {}
        self.flow = self.generate_flow(**flow_state)
        self.ctx['get_flow'].return_value = self.flow
        self.initial_ticks = 1
        self.initial_state = {
            'data': {'ticks': self.initial_ticks,
                     'flow_uuid': self.flow['uuid']},
            'input': {'flow_type': 'some flow type'},
            'status': 'RUNNING'}
        self.task = self.generate_task(**self.initial_state)
        self.task.tick(ctx=self.ctx)

class IncompleteFlowTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(flow_state={'status': 'PENDING'})

    def test_has_expected_state(self):
        self.assertEqual(self.task.status, self.initial_state['status'])
        self.assertEqual(self.task.data, {**self.initial_state['data'],
                                          'ticks': self.initial_ticks + 1})

class CompletedFlowTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(flow_state={
            'status': 'COMPLETED',
            'data': 'some flow_data',
            'output': 'some output',
        })

    def test_has_expected_state(self):
        self.assertEqual(self.task.status, 'COMPLETED')
        self.assertEqual(self.task.data, {**self.initial_state['data'],
                                          'ticks': self.initial_ticks + 1})
        self.assertEqual(self.task.output, self.flow['output'])

class FailedFlowTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(flow_state={
            'status': 'FAILED',
            'error': 'some error'
        })

    def test_has_expected_state(self):
        self.assertEqual(self.task.status, 'FAILED')
        self.assertEqual(self.task.data, {**self.initial_state['data'],
                                          'ticks': self.initial_ticks + 1})
        self.assertEqual(self.task.error, self.flow['error'])

if __name__ == '__main__':
    unittest.main()
