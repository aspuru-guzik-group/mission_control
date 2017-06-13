from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import flow_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = MagicMock()
        self.task_context = self.generate_task_context()
        self.task_handler = flow_task_handler.FlowTaskHandler()

    def generate_flow(self, uuid=None, status='PENDING', **flow_state):
        flow = defaultdict(MagicMock)
        flow.update({'status': status, **flow_state})
        return flow

    def generate_task_context(self, **kwargs):
        task_context = {'task': self.task, **kwargs}
        return task_context

    def setup_mocks(self, attrs=None):
        patchers = {attr: patch.object(flow_task_handler, attr)
                    for attr in attrs}
        mocks = {}
        for key, patcher in patchers.items():
            self.addCleanup(patcher.stop)
            mocks[key] = patcher.start()
        return mocks

    def setup_task_handler_mocks(self, attrs=None):
        for attr in attrs: setattr(self.task_handler, attr, MagicMock())

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.initial_task = defaultdict(MagicMock, **{'data': {}})
        self.task = {**self.initial_task}
        self.setup_task_handler_mocks(attrs=['create_flow'])
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_initial_tick_creates_flow(self):
        self.assertEqual(self.task_handler.create_flow.call_args,
                         call(task=self.task, task_context=self.task_context))

    def test_has_flow_meta(self):
        self.assertEqual(self.task['data']['flow_meta'],
                         self.task_handler.create_flow.return_value)

class CreateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mocks = self.setup_mocks(attrs=['flow_engine'])
        self.result = self.task_handler.create_flow(
            task=self.task, task_context=self.task_context)

    def test_creates_flow(self):
        self.assertEqual(
            self.mocks['flow_engine'].FlowEngine.generate_flow.call_args,
            call(flow_spec=self.task['task_params']['flow_spec'])
        )

    def test_creates_flow_serialization(self):
        self.assertEqual(
            self.mocks['flow_engine'].FlowEngine.serialize_flow.call_args,
            call(flow=(self.mocks['flow_engine'].FlowEngine.generate_flow
                       .return_value))
        )

    def test_returns_flow_meta(self):
        expected_result = {
            'serialization': (self.mocks['flow_engine'].FlowEngine
                              .serialize_flow.return_value)
        }
        self.assertEqual(self.result, expected_result)

class IntermediateTickMixin(object):
    def do_intermediate_tick(self, flow_state=None):
        if not flow_state: flow_state = {}
        self.flow = MagicMock()
        for attr, value in flow_state.items(): setattr(self.flow, attr, value)
        self.task_handler.get_flow = MagicMock(return_value=self.flow)
        self.initial_task = {
            'data': {'flow_uuid': self.flow.uuid},
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

    def test_has_expected_data(self):
        self.assertEqual(self.task['data']['flow_data'], self.flow.data)

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
