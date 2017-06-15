from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock
import uuid

from ..job_task_handler import JobTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = defaultdict(MagicMock, **{'data': {}})
        self.job_context = {
            'mc.tasks.job.create_job': MagicMock(),
            'mc.tasks.job.get_job': MagicMock()
        }
        self.task_ctx = self.generate_task_ctx()
        self.task_handler = JobTaskHandler(task_ctx=self.task_ctx)

    def generate_job(self, key=None, status='PENDING', **job_state):
        job = defaultdict(MagicMock)
        if not key: key = str(uuid.uuid4())
        job.update({'key': key, 'status': status, **job_state})
        return job

    def generate_task_ctx(self, **kwargs):
        return {'task': self.task, **self.job_context, **kwargs}

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.create_job = MagicMock()
        self.task_handler.initial_tick()

    def test_creates_job(self):
        self.assertEqual(self.task_handler.create_job.call_args, call())

    def test_task_stores_job_meta(self):
        self.assertEqual(self.task_handler.task['data']['_job_task_job_meta'],
                         self.task_handler.create_job.return_value)

class CreateJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.result = self.task_handler.create_job()

    def test_dispatches_to_ctx_fn(self):
        self.assertEqual(
            self.task_handler.task_ctx['mc.tasks.job.create_job'].call_args,
            call(job_kwargs={
                'job_spec': self.task['task_params'].get('job_spec')
            })
        )
        self.assertEqual(
            self.result,
            self.task_ctx['mc.tasks.job.create_job'].return_value
        )

class IntermediateTickMixin(object):
    def do_intermediate_tick(self, job_state=None):
        if not job_state: job_state = {}
        self.job = defaultdict(MagicMock)
        for k, v in job_state.items(): self.job[k] = v
        self.task_handler.get_job = MagicMock(return_value=self.job)
        self.initial_task = {
            'data': {'_job_task_job_meta': self.job['key']},
            'task_params': {'job_spec': 'some job spec'},
            'status': 'some_status'
        }
        self.task_handler.task = defaultdict(MagicMock, **self.initial_task)
        self.task_handler.intermediate_tick()

class IncompleteJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.do_intermediate_tick(job_state={'status': 'PENDING'})

    def test_has_expected_state(self):
        self.assertEqual(self.task_handler.task['status'],
                         self.initial_task['status'])
        self.assertEqual(self.task_handler.task['data'],
                         self.initial_task['data'])

class CompletedJobTestCase(BaseTestCase, IntermediateTickMixin):
    def setUp(self):
        super().setUp()
        self.do_intermediate_tick(job_state={'status': 'COMPLETED',
                                             'data': MagicMock()})

    def test_has_expected_status(self):
        self.assertEqual(self.task_handler.task['status'], 'COMPLETED')

    def test_has_expected_artifact(self):
        self.assertEqual(self.task_handler.task['data']['artifact'],
                         self.job.get('data').get('artifact'))

    def test_has_expected_stdout(self):
        self.assertEqual(self.task_handler.task['data']['std_logs'],
                         self.job.get('data').get('std_logs'))

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
