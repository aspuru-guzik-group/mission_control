import unittest
from unittest.mock import call, MagicMock

from ..nodes.flow_wrapper import FlowWrapperNode


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()

    def generate_node(self, **node_kwargs):
        return FlowWrapperNode(**node_kwargs)

class InitializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_state = {
            'status': 'some status',
            'data': {'some': 'data'}
        }
        self.node_data = {'_wrapped_flow_state': self.flow_state}
        self.node = self.generate_node(data=self.node_data)

    def test_deserializes_wrapped_flow_from_node_data(self):
        self.assertEqual(self.node.wrapped_flow.status,
                         self.flow_state['status'])
        self.assertEqual(self.node.wrapped_flow.data,
                         self.flow_state['data'])

class IntermediateTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()
        self.node.data['ticks'] = 1
        self.node.wrapped_flow = MagicMock()
        self.node.on_flow_completed = MagicMock()
        self.node.on_flow_failed = MagicMock()

    def test_ticks_wrapped_flow(self):
        ctx = MagicMock()
        self.node.tick(ctx=ctx)
        self.assertEqual(self.node.wrapped_flow.tick.call_args,
                         call(ctx=ctx))

    def test_sets_wrapped_flow_state(self):
        self.node.tick()
        expected_wrapped_flow_state = {
            'status': self.node.wrapped_flow.status,
            'data': self.node.wrapped_flow.data,
        }
        self.assertEqual(self.node.data['_wrapped_flow_state'],
                         expected_wrapped_flow_state)

    def test_calls_on_wrapped_flow_completed(self):
        self.node.wrapped_flow.status = 'COMPLETED'
        self.node.tick()
        self.assertEqual(self.node.on_wrapped_flow_completed.call_count, 1)

    def test_calls_on_wrapped_flow_failed(self):
        self.node.wrapped_flow.status = 'FAILED'
        self.node.tick()
        self.assertEqual(self.node.on_wrapped_flow_failed.call_count, 1)

    def test_propagates_status_otherwise(self):
        self.node.tick()
        self.assertEqual(self.node.status, self.node.wrapped_flow.status)

class OnWrappedFlowCompletedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()
        self.node.wrapped_flow = MagicMock()
        self.node.job_node.data = {'output': 'some output'}
        self.node.on_wrapped_flow_completed()

    def test_copies_wrapped_flow_output_to_node_output(self):
        self.assertEqual(
            self.node.data['output'], self.node.wrapped_flow.data['output'])

    def test_sets_status(self):
        self.assertEqual(self.node.status, 'COMPLETED')

class OnWrappedFlowFailedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()
        self.node.wrapped_flow = MagicMock()
        self.node.wrapped_flow.data = {'error': 'some error'}
        self.node.on_wrapped_flow_failed()

    def test_copies_wrapped_flow_error_to_node_error(self):
        self.assertEqual(
            self.node.data['error'], self.node.wrapped_flow.data['error'])

    def test_sets_status(self):
        self.assertEqual(self.node.status, 'FAILED')

if __name__ == '__main__':
    unittest.main()
