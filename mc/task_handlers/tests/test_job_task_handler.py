from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock
import uuid

from ..job_task_handler import JobTaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task = MagicMock()
        self.job_context = {
            'mc.tasks.job.create_job': MagicMock(),
            'mc.tasks.job.get_job': MagicMock()
        }
        self.task_handler = JobTaskHandler()

    def generate_job(self, key=None, status='PENDING', **job_state):
        job = defaultdict(MagicMock)
        if not key: key = str(uuid.uuid4())
        job.update({'key': key, 'status': status, **job_state})
        return job

    def generate_task_context(self, **kwargs):
        return {'task': self.task, **self.job_context, **kwargs}

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.initial_task = {'data': {}}
        self.task = {**self.initial_task}
        self.task_context = self.generate_task_context()
        self.task_handler.create_job = MagicMock()
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_creates_job(self):
        self.assertEqual(self.task_handler.create_job.call_args,
                         call(task=self.task, task_context=self.task_context))

    def test_task_stores_job_meta(self):
        self.assertEqual(self.task['data']['_job_meta'],
                         self.task_handler.create_job.return_value)

class CreateJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task = {
            'data': {},
            'task_params': {'job_spec': 'some job_spec'}
        }
        self.task_context = self.generate_task_context()
        self.result = self.task_handler.create_job(
            task=self.task, task_context=self.task_context)

    def test_dispatches_to_ctx_fn(self):
        self.assertEqual(
            self.task_context['mc.tasks.job.create_job'].call_args,
            call(job_kwargs={'job_spec': self.task['task_params']['job_spec']})
        )
        self.assertEqual(
            self.result,
            self.task_context['mc.tasks.job.create_job'].return_value
        )

class IntermediateTickMixin(object):
    def do_intermediate_tick(self, job_state=None):
        if not job_state: job_state = {}
        self.job = defaultdict(MagicMock)
        for k, v in job_state.items(): self.job[k] = v
        self.task_handler.get_job = MagicMock(return_value = self.job)
        self.initial_task = {
            'data': {'_job_meta': self.job['key']},
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

    def test_has_expected_artifact(self):
        self.assertEqual(self.task['data']['artifact'],
                         self.job['data'].get('artifact'))

    def test_has_expected_stdout(self):
        self.assertEqual(self.task['data']['stdout'],
                         self.job['data'].get('stdout'))

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
