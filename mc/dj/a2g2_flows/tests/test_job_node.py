import unittest
from unittest.mock import call, MagicMock
from uuid import uuid4
from django.test import TestCase

from ..nodes.job import JobNode


class BaseTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.ctx = {
            'create_job': MagicMock()
        }

    def generate_node(self, **node_kwargs):
        return JobNode(**node_kwargs)

    def generate_job(self, id=None, status='PENDING', **job_kwargs):
        if not id: id = str(uuid4())
        return {'id': id, 'status': status, **job_kwargs}

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.initial_data = {
            'input': {
                'job_type': 'some job_type',
                'job_spec': 'some job_spec',
            }
        }
        self.node = self.generate_node(data=self.initial_data)
        self.node.tick(ctx=self.ctx)

    def test_initial_tick_creates_job(self):
        expected_call_args = call(job_kwargs={
            'type': self.initial_data['input']['job_type'],
            'spec': self.initial_data['input']['job_spec'],
        })
        self.assertEqual(self.ctx['create_job'].call_args, expected_call_args)

    def test_has_expected_data(self):
        expected_data = {**self.initial_data, 
                         'job_id': self.ctx['create_job'].return_value,
                         'ticks': 1}
        self.assertEqual(self.node.data, expected_data)

    def test_has_expected_status(self):
        self.assertEqual(self.node.status, 'RUNNING')

class IntermediateTickMixin(object):
    def setup_for_intermediate_tick(self, job_kwargs=None):
        if not job_kwargs: job_kwargs = {}
        self.job = self.generate_job(**job_kwargs)
        self.initial_ticks = 1
        self.initial_data = {'ticks': self.initial_ticks,
                              'job_id': self.job['id'],
                              'status': 'RUNNING'}
        self.node = self.generate_node(data=self.initial_data,
                                       jobs={self.job['id']: self.job})
        self.node.tick()

class IncompleteJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(job_kwargs={'status': 'PENDING'})


    def test_has_expected_data(self):
        self.assertEqual(self.node.status, 'RUNNING')
        expected_data = {
            **self.initial_data,
            'ticks': self.initial_ticks + 1,
        }
        self.assertEqual(self.node.data, expected_data)

class CompletedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(job_kwargs={
            'status': 'COMPLETED',
            'data': 'some job_data',
        })

    def test_has_expected_status(self):
        self.assertEqual(self.node.status, 'COMPLETED')

    def test_has_expected_data(self):
        expected_data = {
            **self.initial_data,
            'ticks': self.initial_ticks + 1,
            'outputs': {
                'job_data': self.job['data'],
            }
        }
        self.assertEqual(self.node.data, expected_data)

class FailedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(job_kwargs={
            'status': 'FAILED',
            'data': {'error': 'some error'},
        })

    def test_has_expected_data(self):
        self.assertEqual(self.node.status, 'FAILED')
        expected_data = {
            **self.initial_data,
            'ticks': self.initial_ticks + 1,
            'error': self.job['data']['error'],
        }
        self.assertEqual(self.node.data, expected_data)

if __name__ == '__main__':
    unittest.main()
