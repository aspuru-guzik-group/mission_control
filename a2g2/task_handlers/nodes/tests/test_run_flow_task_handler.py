from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock
from uuid import uuid4

from .. import run_flow_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = MagicMock()
        self.flow_ctx = {
            'create_flow': MagicMock(),
            'get_flow': MagicMock()
        }
        self.task_context = self.generate_task_context()
        self.task_handler = run_flow_task_handler.RunFlowTaskHandler()

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
        self.initial_task = defaultdict(MagicMock, **{'data': {}})
        self.task = {**self.initial_task}
        self.task_handler.create_flow = MagicMock()
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_initial_tick_creates_flow(self):
        self.assertEqual(self.task_handler.create_flow.call_args,
                         call(task=self.task, task_context=self.task_context))

    def test_has_flow_uuid(self):
        self.assertEqual(self.task['data']['flow_uuid'],
                         self.task_handler.create_flow.return_value['uuid'])

class CreateFlowTestCase(BaseTestCase):
    def test_creates_flow_w_expected_kwargs(self):
        self.task_handler.generate_flow_kwargs = MagicMock()
        flow = self.task_handler.create_flow(task=self.task,
                                             task_context=self.task_context)
        expected_flow_kwargs = \
                self.task_handler.generate_flow_kwargs.return_value
        self.assertEqual(self.flow_ctx['create_flow'].call_args,
                         call(flow_kwargs=expected_flow_kwargs))
        self.assertEqual(self.task_handler.generate_flow_kwargs.call_args,
                         call(flow_spec=self.task['task_params']['flow_spec']))
        self.assertEqual(flow, self.flow_ctx['create_flow'].return_value)

class GenerateFlowKwargsTestCase(BaseTestCase):
    @patch.object(run_flow_task_handler, 'flow_engine')
    def test_generate_expected_flow_kwargs(self, mock_flow_engine):
        expected_FlowEngine = mock_flow_engine.FlowEngine
        flow_spec = MagicMock()
        flow_kwargs = self.task_handler.generate_flow_kwargs(
            flow_spec=flow_spec)
        expected_flow_kwargs = {
            'label': flow_spec.get('label'),
            'serialization': expected_FlowEngine.serialize_flow.return_value
        }
        self.assertEqual(flow_kwargs, expected_flow_kwargs)
        self.assertEqual(
            expected_FlowEngine.serialize_flow.call_args,
            call(flow=expected_FlowEngine.generate_flow.return_value))
        self.assertEqual(expected_FlowEngine.generate_flow.call_args,
                         call(flow_spec=flow_spec))

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
