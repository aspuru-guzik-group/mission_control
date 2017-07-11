import unittest
from unittest.mock import call, MagicMock, patch

from .. import tick


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.subcommand = tick.TickSubcommand(utils=MagicMock())

    def mockify_subcommand_attrs(self, attrs=None):
        for attr in attrs: setattr(self.subcommand, attr, MagicMock())

    def mockify_module_attrs(self, attrs=None, module=tick):
        module_mocks = {}
        for attr in attrs: 
            patcher = patch.object(module, attr)
            self.addCleanup(patcher.stop)
            module_mocks[attr] = patcher.start()
        return module_mocks

class _RunTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_tick_flow_runner',
                                             '_tick_job_runner',
                                             '_tick_jobman'])
        self.subcommand._run()

    def test_ensures_db(self):
        self.assertEqual(self.subcommand.utils.ensure_db.call_args, call())

    def test_ensures_queues(self):
        self.assertEqual(self.subcommand.utils.ensure_queues.call_args, call())

    def test_ticks_flow_runner(self):
        self.assertEqual(self.subcommand._tick_flow_runner.call_args, call())

    def test_ticks_job_runner(self):
        self.assertEqual(self.subcommand._tick_job_runner.call_args, call())

    def test_ticks_jobman(self):
        self.assertEqual(self.subcommand._tick_jobman.call_args, call())

class _TickFlowRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.subcommand._tick_flow_runner()

    def test_ticks_flow_runner(self):
        self.assertEqual(self.subcommand.utils.flow_runner.tick.call_args,
                         call())

class _TickJobRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.subcommand._tick_job_runner()

    def test_ticks_job_runner(self):
        self.assertEqual(self.subcommand.utils.job_runner.tick.call_args,
                         call())

class _TickJobManTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.subcommand._tick_jobman()

    def test_ticks_jobman(self):
        self.assertEqual(self.subcommand.utils.jobman.tick.call_args, call())
