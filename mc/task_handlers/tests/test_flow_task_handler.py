from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from .. import flow_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = defaultdict(MagicMock, {'data': defaultdict(MagicMock)})
        self.task_ctx = self.generate_task_ctx()
        self.task_handler = flow_task_handler.FlowTaskHandler(
            task_ctx=self.task_ctx)

    def generate_task_ctx(self, **kwargs):
        return defaultdict(MagicMock, **{'task': self.task, **kwargs})

    def setup_handler_mocks(self, attrs=None):
        for attr in attrs: setattr(self.task_handler, attr, MagicMock())

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_handler_mocks(attrs=['create_flow_record'])
        self.task_handler.initial_tick()

    def test_initial_tick_creates_flow_record(self):
        self.assertEqual(self.task_handler.create_flow_record.call_args, call())

    def test_stores_flow_meta(self):
        self.assertEqual(self.task['data']['_flow_task_flow_meta'],
                         self.task_handler.create_flow_record.return_value)

class CreateFlowRecordTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_handler_mocks(attrs=['create_flow', 'flow_record_client',
                                        'flow_engine'])
        self.result = self.task_handler.create_flow_record()

    def test_creates_flow(self):
        self.assertEqual(self.task_handler.create_flow.call_args, call())

    def test_calls_flow_to_flow_dict(self):
        self.assertEqual(
            self.task_handler.flow_engine.flow_to_flow_dict.call_args,
            call(flow=self.task_handler.create_flow.return_value)
        )

    def test_calls_create_flow_record(self):
        self.assertEqual(
            self.task_handler.flow_record_client.create_flow_record.call_args,
            call(flow_kwargs=(self.task_handler.flow_engine.flow_to_flow_dict
                              .return_value))
        )

    def test_returns_flow_meta(self):
        self.assertEqual(
            self.result,
            {'key': (self.task_handler.flow_record_client.create_flow_record
                     .return_value['key'])}
        )

class CreateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.expected_flow = (self.task_handler.flow_engine.flow_spec_to_flow
                              .return_value)
        self.expected_flow.data = {}
        self.result = self.task_handler.create_flow()

    def test_calls_flow_engine_flow_spec_to_flow(self):
        self.assertEqual(
            self.task_handler.flow_engine.flow_spec_to_flow.call_args,
            call(flow_spec=self.task['task_params']['flow_spec'])
        )

    def test_sets_depth(self):
        self.assertEqual(self.expected_flow.depth,
                         self.task_ctx['flow'].depth + 1)

    def test_sets_parent_key(self):
        self.assertEqual(self.expected_flow.data['parent_key'],
                         self.task_ctx['flow'].key)

    def test_returns_flow(self):
        self.assertEqual(self.result, self.expected_flow)

class IntermediateTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_handler_mocks(attrs=['get_flow', 'handle_flow_status'])
        self.task_handler.intermediate_tick()

    def test_gets_flow(self):
        self.assertEqual(self.task_handler.get_flow.call_args, call())

    def test_handles_flow_status(self):
        self.assertEqual(self.task_handler.handle_flow_status.call_args,
                         call(flow=self.task_handler.get_flow.return_value))

class GetFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_handler_mocks(attrs=['get_flow_record'])
        self.result = self.task_handler.get_flow()

    def test_calls_get_flow_record(self):
        self.assertEqual(self.task_handler.get_flow_record.call_args, call())

    def test_returns_engine_flow_dict_to_flow(self):
        self.assertEqual(
            self.task_handler.flow_engine.flow_dict_to_flow.call_args,
            call(flow_dict=self.task_handler.get_flow_record.return_value)
        )
        self.assertEqual(
            self.result,
            self.task_handler.flow_engine.flow_dict_to_flow.return_value
        )

class GetFlowRecordTestCase(BaseTestCase):
    def test_dispatches_to_flow_record_client(self):
        result = self.task_handler.get_flow_record()
        self.assertEqual(
            self.task_handler.flow_record_client.get_flow_record.call_args,
            call(flow_meta=self.task['data']['_flow_task_flow_meta'])
        )
        self.assertEqual(
            result,
            self.task_handler.flow_record_client.get_flow_record.return_value
        )
