from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock

from .. import inline_flow_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = defaultdict(MagicMock, **{'data': {}})
        self.flow_engine = MagicMock()
        self.task_ctx = self.generate_task_ctx(flow_engine=self.flow_engine)
        self.task_handler = inline_flow_task_handler.InlineFlowTaskHandler(
            task_ctx=self.task_ctx)

    def generate_flow(self, uuid=None, status='PENDING', **flow_state):
        flow = defaultdict(MagicMock)
        flow.update({'status': status, **flow_state})
        return flow

    def generate_task_ctx(self, **kwargs):
        task_ctx = defaultdict(MagicMock, **{'task': self.task, **kwargs})
        return task_ctx

    def setup_mocks(self, attrs=None):
        patchers = {attr: patch.object(inline_flow_task_handler, attr)
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
        self.setup_task_handler_mocks(attrs=['create_flow',
                                             'tick_flow_until_has_no_pending'])
        self.task_handler.initial_tick()

    def test_creates_flow(self):
        self.assertEqual(self.task_handler.create_flow.call_args(), call())

    def test_ticks_created_flow(self):
        self.assertEqual(
            self.task_handler.tick_flow_until_has_no_pending.call_args,
            call(flow=self.task_handler.create_flow.return_value))

class CreateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.result = self.task_handler.create_flow()

    def test_dispatches_to_flow_engine(self):
        self.assertEqual(
            self.flow_engine.generate_flow.call_args,
            call(flow_spec=self.task_handler.task['task_params']['flow_spec']))
        self.assertEqual(self.result,
                         self.flow_engine.generate_flow.return_value)

class TickFlowUntilHasNoPendingTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = MagicMock()
        self.setup_task_handler_mocks(attrs=['persist_flow',
                                             'handle_completed_flow'])

    def _tick(self):
        self.task_handler.tick_flow_until_has_no_pending(flow=self.flow)

    def test_dispatches_to_flow_engine(self):
        self._tick()
        expected_task_ctx = {
            'parent_task_ctx': self.task_ctx,
            **{k: v for k, v in self.task_ctx.items()
               if k not in {'flow', 'task'}}
        }
        self.assertEqual(
            self.flow_engine.tick_flow_until_has_no_pending.call_args,
            call(flow=self.flow, task_ctx=expected_task_ctx))

    def test_persists_flow(self):
        self._tick()
        self.assertEqual(self.task_handler.persist_flow.call_args,
                         call(flow=self.flow))

    def test_handles_completed_flow(self):
        def mock_tick_flow(*args, flow=None, **kargs): flow.status = 'COMPLETED'
        self.flow_engine.tick_flow_until_has_no_pending.side_effect = \
                mock_tick_flow
        self._tick()
        self.assertEqual(self.task_handler.handle_completed_flow.call_args,
                         call(flow=self.flow))

class PersistFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = MagicMock()
        self.task_handler.persist_flow(flow=self.flow)

    def test_serializes_flow(self):
        self.assertEqual(self.flow_engine.serialize_flow.call_args,
                         call(flow=self.flow))

    def test_saves_flow_meta_to_task_data(self):
        self.assertEqual(self.task_handler.task['data']['_flow_task_flow_meta'],
                         {'serialization': (self.flow_engine.serialize_flow
                                            .return_value)})
            
    def test_sets_task_status(self):
        self.assertEqual(self.task_handler.task['status'], 'RUNNING')

class HandleCompletedFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow = MagicMock()

    def _handle(self):
        self.task_handler.handle_completed_flow(flow=self.flow)

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
