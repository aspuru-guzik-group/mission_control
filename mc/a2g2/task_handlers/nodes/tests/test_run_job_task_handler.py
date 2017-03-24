import unittest
from unittest.mock import call, MagicMock
from uuid import uuid4

from ..run_job_task_handler import RunJobTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = MagicMock()
        self.flow_ctx = {
            'create_job': MagicMock(),
            'get_job': MagicMock()
        }
        self.task_handler = RunJobTaskHandler()

    def generate_job(self, uuid=None, status='PENDING', **job_state):
        if not uuid: uuid = str(uuid4())
        return {'uuid': uuid, 'status': status, **job_state}

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
        self.initial_task = {
            'data': {},
            'task_params': {
                'job_spec': 'some job_spec'
            }
        }
        self.task = {**self.initial_task}
        self.task_handler.initial_tick(
            task=self.task,
            task_context=self.generate_task_context()
        )

    def test_initial_tick_creates_job(self):
        self.assertEqual(
            self.flow_ctx['create_job'].call_args,
            call(job_kwargs={
                'job_spec': self.initial_task['task_params']['job_spec'],
            })
        )

    def test_has_job_uuid(self):
        self.assertEqual(self.task['data']['job_uuid'],
                         self.flow_ctx['create_job'].return_value['uuid'])

class IntermediateTickMixin(object):
    def do_intermediate_tick(self, job_state=None):
        if not job_state: job_state = {}
        self.job = self.generate_job(**job_state)
        self.flow_ctx['get_job'].return_value = self.job
        self.initial_task = {
            'data': {'job_uuid': self.job['uuid']},
            'task_params': {'job_spec': 'some job spec'},
            'status': 'some_status'
        }
        self.task = {**self.initial_task}
        self.task_handler.intermediate_tick(
            task=self.task,
            task_context=self.generate_task_context()
        )

class IncompleteJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.do_intermediate_tick(job_state={'status': 'PENDING'})

    def test_has_expected_state(self):
        self.assertEqual(self.task['status'], self.initial_task['status'])
        self.assertEqual(self.task['data'], self.initial_task['data'])

class CompletedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.do_intermediate_tick(job_state={'status': 'COMPLETED'})

    def test_has_expected_status(self):
        self.assertEqual(self.task['status'], 'COMPLETED')

class FailedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def test_throws_exception(self):
        with self.assertRaises(Exception) as ctx:
            error = 'some error'
            self.do_intermediate_tick(job_state={
                'status': 'FAILED',
                'data': {'error': error},
            })
            self.assertTrue(error in ctx.exception)

if __name__ == '__main__':
    unittest.main()
