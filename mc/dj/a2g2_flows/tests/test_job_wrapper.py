import unittest
from unittest.mock import call, MagicMock

from ..nodes.job_wrapper import JobWrapperNode


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        class TestNode(JobWrapperNode):
            def get_job_input(self):  return {'job_type': 'some job type'}
        self.TestNode = TestNode

    def generate_node(self, **node_kwargs):
        return self.TestNode(**node_kwargs)

class InitializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.nested_job_state = {
            'status': 'some status',
            'data': {'some': 'data'}
        }
        self.wrapper_node_data = {'_job_node_state': self.nested_job_state}
        self.node = self.generate_node(data=self.wrapper_node_data)

    def test_job_node_gets_state_from_nested_job_state(self):
        self.assertEqual(
            self.node.job_node.status, self.nested_job_state['status'])
        self.assertEqual(self.node.job_node.data, self.nested_job_state['data'])

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()

    def test_job_node_gets_input_data(self):
        self.node.tick()
        self.assertEqual(
            self.node.job_node.data['input'], self.node.get_job_input())

class IntermediateTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()
        self.node.data['ticks'] = 1
        self.node.job_node = MagicMock()
        self.node.on_job_completed = MagicMock()
        self.node.on_job_failed = MagicMock()

    def test_calls_job_node_tick(self):
        ctx = MagicMock()
        self.node.tick(ctx=ctx)
        self.assertEqual(self.node.job_node.tick.call_args, call(ctx=ctx))

    def test_sets_job_node_state(self):
        self.node.tick()
        expected_job_node_state = {'status': self.node.job_node.status,
                                   'data': self.node.job_node.data}
        self.assertEqual(
            self.node.data['_job_node_state'], expected_job_node_state)

    def test_calls_on_job_completed(self):
        self.node.job_node.status = 'COMPLETED'
        self.node.tick()
        self.assertEqual(self.node.on_job_completed.call_count, 1)

    def test_calls_on_job_failed(self):
        self.node.job_node.status = 'FAILED'
        self.node.tick()
        self.assertEqual(self.node.on_job_failed.call_count, 1)

    def test_propagates_status_otherwise(self):
        self.node.tick()
        self.assertEqual(self.node.status, self.node.job_node.status)

class OnJobCompletedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()
        self.node.job_node = MagicMock()
        self.node.job_node.data = {'output': 'some output'}
        self.node.on_job_completed()

    def test_copies_job_output_to_node_output(self):
        self.assertEqual(
            self.node.data['output'], self.node.job_node.data['output'])

    def test_sets_status(self):
        self.assertEqual(self.node.status, 'COMPLETED')

class OnJobFailedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node = self.generate_node()
        self.node.job_node = MagicMock()
        self.node.job_node.data = {'error': 'some error'}
        self.node.on_job_failed()

    def test_copies_job_error_to_node_error(self):
        self.assertEqual(
            self.node.data['error'], self.node.job_node.data['error'])

    def test_sets_status(self):
        self.assertEqual(self.node.status, 'FAILED')

if __name__ == '__main__':
    unittest.main()
