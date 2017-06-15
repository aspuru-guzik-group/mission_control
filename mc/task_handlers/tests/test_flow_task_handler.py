from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from .. import flow_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = defaultdict(MagicMock, {'data': defaultdict(MagicMock)})
        self.flow_engine = MagicMock()
        self.task_ctx = self.generate_task_ctx()
        self.task_handler = flow_task_handler.FlowTaskHandler(
            task_ctx=self.task_ctx)

    def generate_task_ctx(self, **kwargs):
        task_ctx = defaultdict(
            MagicMock,
            **{'task': self.task, 'flow_engine': self.flow_engine, **kwargs}
        )
        return task_ctx

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.create_flow = MagicMock()
        self.task_handler.initial_tick()

    def test_initial_tick_creates_flow(self):
        self.assertEqual(self.task_handler.create_flow.call_args, call())

    def test_has_flow_meta(self):
        self.assertEqual(self.task['data']['_flow_task_flow_meta'],
                         self.task_handler.create_flow.return_value)

class CreateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.generate_flow_record_kwargs = MagicMock()
        self.result = self.task_handler.create_flow()

    def test_flow_record_generates_kwargs(self):
        self.assertEqual(
            self.task_handler.generate_flow_record_kwargs.call_args, call())

    def test_calls_create_flow_record_ctx_fn(self):
        self.assertEqual(
            self.task_ctx['mc.tasks.flow.create_flow_record'].call_args,
            call(flow_kwargs=(self.task_handler.generate_flow_record_kwargs
                              .return_value))
        )

    def test_returns_flow_record_meta(self):
        self.assertEqual(
            self.result,
            {'key': (self.task_ctx['mc.tasks.flow.create_flow_record']
                     .return_value)['key']}
        )

class GenerateFlowKwargsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.serialize_flow_spec = MagicMock()
        self.result = self.task_handler.generate_flow_record_kwargs()

    def test_gets_serialization_for_spec(self):
        self.assertEqual(
            self.task_handler.serialize_flow_spec.call_args,
            call(flow_spec=self.task['task_params']['flow_spec']))

    def test_returns_expected_kwargs(self):
        expected_kwargs = {
            'label': self.task['task_params']['flow_spec'].get('label'),
            'serialization': self.task_handler.serialize_flow_spec.return_value
        }
        self.assertEqual(self.result, expected_kwargs)

class GetSerializationForFlowSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_spec = MagicMock()
        self.result = self.task_handler.serialize_flow_spec(
            flow_spec=self.flow_spec)

    def test_gets_serialization(self):
        self.assertEqual(self.flow_engine.generate_flow.call_args,
                         call(flow_spec=self.flow_spec))
        self.assertEqual(self.flow_engine.serialize_flow.call_args,
                         call(flow=self.flow_engine.generate_flow.return_value))
        self.assertEqual(self.result,
                         self.flow_engine.serialize_flow.return_value)

class IntermediateTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.get_flow = MagicMock()
        self.task_handler.handle_flow_status = MagicMock()
        self.task_handler.intermediate_tick()

    def test_gets_flow(self):
        self.assertEqual(self.task_handler.get_flow.call_args, call())

    def test_handles_flow_status(self):
        self.assertEqual(self.task_handler.handle_flow_status.call_args,
                         call(flow=self.task_handler.get_flow.return_value))

class GetFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.deserialize_flow_record = MagicMock()
        self.result = self.task_handler.get_flow()

    def test_calls_get_flow_ctx_fn(self):
        self.assertEqual(
            self.task_ctx['mc.tasks.flow.get_flow_record'].call_args,
            call(flow_meta=self.task['data']['_flow_task_flow_meta'])
        )

    def test_deserializes_ctx_fn_result(self):
        self.assertEqual(
            self.task_handler.deserialize_flow_record.call_args,
            call(flow_record=(self.task_ctx['mc.tasks.flow.get_flow_record']
                              .return_value))
        )
        self.assertEqual(self.result,
                         self.task_handler.deserialize_flow_record.return_value)

class DeserializeFlowRecordTestCase(BaseTestCase):
    def test_dispatches_to_flow_engine(self):
        flow_record = MagicMock()
        result = self.task_handler.deserialize_flow_record(
            flow_record=flow_record)
        self.assertEqual(self.flow_engine.deserialize_flow.call_args,
                         call(serialized_flow=flow_record['serialization']))
        self.assertEqual(result, self.flow_engine.deserialize_flow.return_value)
