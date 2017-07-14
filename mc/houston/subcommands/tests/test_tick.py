import unittest
from unittest.mock import call, MagicMock, patch

from .. import tick


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.subcommand = tick.TickSubcommand(utils=MagicMock())
        self.mockify_subcommand_attrs(attrs=['_sleep'])

    def mockify_subcommand_attrs(self, attrs=None):
        for attr in attrs: setattr(self.subcommand, attr, MagicMock())

    def mockify_module_attrs(self, attrs=None, module=tick):
        if not hasattr(self, 'modmocks'): self.modmocks = {}
        for attr in attrs: 
            patcher = patch.object(module, attr)
            self.addCleanup(patcher.stop)
            self.modmocks[attr] = patcher.start()

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

class NTicksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_tick'])
        self.subcommand.kwargs['nticks'] = 3
        self.subcommand._run()

    def test_ticks_for_nticks(self):
        self.assertEqual(len(self.subcommand._tick.call_args_list),
                         self.subcommand.kwargs['nticks'])

class MaxTicksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_tick'])
        self.subcommand.kwargs['nticks'] = 3
        self.subcommand.kwargs['max_ticks'] = \
                self.subcommand.kwargs['nticks'] - 1

    def test_raises_if_exceeds_max_ticks(self):
        with self.assertRaises(Exception):
            self.subcommand._run()

class IntervalTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_tick'])
        self.subcommand.kwargs['nticks'] = 2
        self.subcommand.kwargs['interval'] = 123
        self.subcommand._run()

    def test_ticks_at_interval(self):
        self.assertEqual(self.subcommand._sleep.call_args,
                         call(self.subcommand.kwargs['interval']))
