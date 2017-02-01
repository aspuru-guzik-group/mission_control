import unittest
from unittest.mock import call, MagicMock
from uuid import uuid4

from ..job_task_engine import JobTaskEngine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.ctx = {
            'create_job': MagicMock(),
            'get_job': MagicMock()
        }
        self.engine = JobTaskEngine()

    def generate_job(self, uuid=None, status='PENDING', **job_state):
        if not uuid: uuid = str(uuid4())
        return {'uuid': uuid, 'status': status, **job_state}

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.initial_task = {
            'input': {
                'job_spec': 'some job_spec'
            }
        }
        self.task = {**self.initial_task}
        self.engine.tick_task(task=self.task, ctx=self.ctx)

    def test_initial_tick_creates_job(self):
        expected_call_args = call(job_kwargs={
            'job_spec': self.initial_task['input']['job_spec'],
        })
        self.assertEqual(self.ctx['create_job'].call_args, expected_call_args)

    def test_has_expected_data(self):
        expected_data = {
            'job_uuid': self.ctx['create_job'].return_value['uuid'],
            'ticks': 1
        }
        self.assertEqual(self.task['data'], expected_data)

    def test_has_expected_status(self):
        self.assertEqual(self.task['status'], 'RUNNING')

class IntermediateTickMixin(object):
    def setup_for_intermediate_tick(self, job_state=None):
        if not job_state: job_state = {}
        self.job = self.generate_job(**job_state)
        self.ctx['get_job'].return_value = self.job
        self.initial_ticks = 1
        self.initial_task = {
            'data': {'ticks': self.initial_ticks,
                     'job_uuid': self.job['uuid']},
            'input': {'job_spec': 'some job spec'},
            'status': 'RUNNING'}
        self.task = {**self.initial_task}
        self.engine.tick_task(task=self.task, ctx=self.ctx)

class IncompleteJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(job_state={'status': 'PENDING'})

    def test_has_expected_state(self):
        self.assertEqual(self.task['status'], self.initial_task['status'])
        self.assertEqual(self.task['data'], {**self.initial_task['data'],
                                             'ticks': self.initial_ticks + 1})

class CompletedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(job_state={
            'status': 'COMPLETED',
            'data': {'output': 'some output'},
        })

    def test_has_expected_state(self):
        self.assertEqual(self.task['status'], 'COMPLETED')
        self.assertEqual(self.task['data'], {**self.initial_task['data'],
                                             'ticks': self.initial_ticks + 1})
        self.assertEqual(self.task['output'], self.job['data']['output'])

class FailedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.setup_for_intermediate_tick(job_state={
            'status': 'FAILED',
            'data': {'error': 'some error'},
        })

    def test_has_expected_state(self):
        self.assertEqual(self.task['status'], 'FAILED')
        self.assertEqual(self.task['data'], {**self.initial_task['data'],
                                             'ticks': self.initial_ticks + 1})
        self.assertEqual(self.task['error'], self.job['data']['error'])

if __name__ == '__main__':
    unittest.main()
