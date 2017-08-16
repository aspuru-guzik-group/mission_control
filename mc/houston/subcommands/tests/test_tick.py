import unittest
from unittest.mock import call, MagicMock, patch

from .. import tick


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.skipTest('RETURN TO THIS LATER!')
        self.subcommand = tick.TickSubcommand(utils=MagicMock())
        self.mockify_subcommand_attrs(attrs=['_sleep'])

    def mockify_subcommand_attrs(self, attrs=None):
        for attr in attrs:
            setattr(self.subcommand, attr, MagicMock())

    def mockify_module_attrs(self, attrs=None, module=tick):
        if not hasattr(self, 'modmocks'):
            self.modmocks = {}
        for attr in attrs:
            patcher = patch.object(module, attr)
            self.addCleanup(patcher.stop)
            self.modmocks[attr] = patcher.start()


class _RunTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_dispatch_on_tickee'])
        self.subcommand._run()

    def test_ensures_db(self):
        self.assertEqual(self.subcommand.utils.ensure_db.call_args, call())

    def test_ensures_queues(self):
        self.assertEqual(self.subcommand.utils.ensure_queues.call_args, call())

    def test_dispatches_on_tickee(self):
        self.assertEqual(self.subcommand._dispatch_on_tickee.call_args, call())


class _DispatchOnTickeeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_get_fn_for_tickee'])
        self.subcommand.kwargs['tickee'] = MagicMock()
        self.result = self.subcommand._dispatch_on_tickee()

    def test_gets_fn_for_tickee(self):
        self.assertEqual(self.subcommand._get_fn_for_tickee.call_args,
                         call(tickee_name=self.subcommand.kwargs['tickee']))

    def test_calls_fn_for_tickee(self):
        self.assertEqual(
            self.subcommand._get_fn_for_tickee.return_value.call_args,
            call()
        )


class _GetFnForTickeeTestCase(BaseTestCase):
    def test_gets_fns_for_tickee_choices(self):
        self.assertEqual(
            self.subcommand.TICKEE_NAMES,
            ['flow', 'flow_runner', 'job_runner', 'jobman', 'all']
        )
        expected_tickee_fns = [
            getattr(self.subcommand, '_run_for_%s_tickee' % tickee_name)
            for tickee_name in self.subcommand.TICKEE_NAMES
        ]
        actual_tickee_fns = [
            self.subcommand._get_fn_for_tickee(tickee_name=tickee_name)
            for tickee_name in self.subcommand.TICKEE_NAMES
        ]
        self.assertEqual(actual_tickee_fns, expected_tickee_fns)


class _RunForFlowTickeeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=[
            '_update_parsed_args_for_flow_tickee', '_run_for_flow_spec',
            '_run_for_flow_key'])

    def _run_for_flow_tickee(self):
        self.subcommand._run_for_flow_tickee()

    def test_updates_parsed_args(self):
        self.assertEqual(
            self.subcommand._update_parsed_args_for_flow_tickee.call_args,
            call()
        )

    def test_runs_for_flow_spec_if_flow_spec(self):
        self.subcommand.kwargs['flow_spec'] = MagicMock()
        self._run_for_flow_tickee()
        self.assertEqual(self.subcommand._run_for_flow_spec.call_args, call())

    def test_runs_for_flow_key_if_flow_key(self):
        self.subcommand.kwargs['flow'] = MagicMock()
        self._run_for_flow_tickee()
        self.assertEqual(self.subcommand._run_for_flow_key.call_args, call())


class _RunForFlowSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_run_for_flow'])
        self.subcommand.kwargs['flow_spec'] = MagicMock()

    def test_converts_flow_to_flow_spec(self):
        self.assertEqual(
            self.subcommand.utils.flow_engine.flow_spec_to_flow.call_args,
            call(flow_spec=self.subcommand.kwargs['flow_spec'])
        )

    def test_runs_for_flow(self):
        self.assertEqual(
            self.subcommand._run_for_flow.call_args,
            call(flow=self.subcommand.flow_engine.flow_spec_to_flow
                 .return_value)
        )


class _RunForFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=[
            '_tick_while', '_get_tick_fn_for_flow',
            '_get_condition_fns_for_flow'])

    def test_gets_tick_fn(self):
        self.assertEqual(self.subcommand._get_tick_fn_for_flow.call_args,
                         call(flow=self.flow))

    def test_gets_condition_fns_for_flow(self):
        call_kwargs = self.subcommand._get_condition_fns_for_flow.call_args[1]
        self.assertTrue(call_kwargs['get_flow'] is not None)

    def test_dispatches_to_tick_while(self):
        self.assertEqual(
            self.subcommand._tick_while.call_args,
            call(
                tick_fn=self.subcommand._get_tick_fn_for_flow.return_value,
                condition_fns=(self.subcommand._get_condition_fns_for_flow
                               .return_value)
            )
        )


class _RunForFlowKeyTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_run_for_flow_runner'])
        self.subcommand.kwargs['flow_key'] = MagicMock()

    def test_calls_run_for_flow_runner_w_flow_key(self):
        self.assertEqual(
            self.subcommand._run_for_flow_runner.call_args,
            call(flow_key=self.subcommand.kwargs['flow_key'])
        )


class _RunForFlowRunnerTickeeTestCase(BaseTestCase):
    def test_something(self):
        self.fail()


class _RunForJobRunnerTickeeTestCase(BaseTestCase):
    def test_something(self):
        self.fail()


class _RunForJobmanTickee(BaseTestCase):
    def test_something(self):
        self.fail()


class _RunForAllTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_subcommand_attrs(attrs=['_tick_flow_runner',
                                             '_tick_job_runner',
                                             '_tick_jobman'])
        self._run_for_all()

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
        self.subcommand.kwargs['max_ticks'] = (
            self.subcommand.kwargs['nticks'] - 1
        )

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
