import unittest
from unittest.mock import call, MagicMock

from ..job_engine import JobEngine


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.default_kwargs = {
            'command_handlers': MagicMock(),
        }
        self.job_engine = self.generate_job_engine()
        self.job = MagicMock()

    def generate_job_engine(self, *args, **kwargs):
        return JobEngine(*args, **{**self.default_kwargs, **kwargs})

class InitTestCase(BaseTestCase):
    def test_sets_command_handlers(self):
        command_handlers = MagicMock()
        job_engine = self.generate_job_engine(command_handlers=command_handlers)
        self.assertEqual(job_engine.command_handlers, command_handlers)

class ExecuteJobTestCase(BaseTestCase):
    def test_dispatches_to_command_handler(self):
        self.job_engine.handle_command = MagicMock()
        self.job_engine.execute_job(job=self.job)
        self.assertEqual(
            self.job_engine.handle_command.call_args,
            call(command=self.job['job_spec']['command'], job=self.job))

class HandleCommandTestCase(BaseTestCase):
    def test_calls_handler_for_command(self):
        self.job_engine.command_handlers = MagicMock()
        command = 'some_command'
        self.job_engine.handle_command(command=command, job=self.job)
        self.assertEqual(self.job_engine.command_handlers[command].call_args,
                         call(job=self.job))
