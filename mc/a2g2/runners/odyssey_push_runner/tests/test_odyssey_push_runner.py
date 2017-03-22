import io
import json
import os
import tempfile
import shlex
import unittest
from unittest.mock import call, DEFAULT, Mock, patch

from django.test import TestCase

from .. import odyssey_push_runner
from ..commands import base as base_command
from ..commands import create_flow_record as create_flow_record_command
from ..commands import run as run_command


class BaseTestCase(TestCase):
    def setUp(self):
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = self.start_patchers()
        self.runner = self.generate_runner()
        self.mock_runner = self.generate_mock_runner()

    def decorate_patchers(self): pass

    def start_patchers(self):
        mocks = {key: patcher.start() for key, patcher in self.patchers.items()}
        return mocks

    def tearDown(self):
        self.stop_patchers()

    def stop_patchers(self):
        for patcher in self.patchers.values(): patcher.stop()

    def assert_is_subdict(self, subdict=None, superdict=None):
        comparison_dict = {key: superdict.get(key, None) for key in subdict}
        self.assertEqual(subdict, comparison_dict)

    def generate_runner(self, *args, run_setup=False, **runner_kwargs):
        return odyssey_push_runner.OdysseyPushRunner(run_setup=run_setup,
                                                     **runner_kwargs)
    def generate_mock_runner(self):
        mock_runner = Mock()
        def call_real_method(method_name, *args, **kwargs):
            return self.call_runner_method_with_mock(
                mock=mock_runner, method_name=method_name, *args, **kwargs)
        mock_runner.call_real_method = call_real_method
        return mock_runner

    def call_runner_method_with_mock(self, mock=None, method_name=None,
                                     *args, **kwargs):
        method =  getattr(odyssey_push_runner.OdysseyPushRunner, method_name)
        return method(mock, *args, **kwargs)

class BaseCommandTestCase(BaseTestCase):
    def setUp(self):
        self.params = {'param1': 'value1', 'param2': 'value2'}
        super().setUp()
        self.setup_mock_add_arguments()

    def decorate_patchers(self):
        self.patchers['BaseCommand'] = patch.multiple(
            base_command.BaseCommand, handle=DEFAULT, add_arguments=DEFAULT)

    def setup_mock_add_arguments(self):
        def mock_add_arguments(parser):
            for param in self.params: parser.add_argument('--%s' % param)
        self.mocks['BaseCommand']['add_arguments'].side_effect = \
                mock_add_arguments

    def params_to_args(self, params=None):
        arg_strs = ["--%s='%s'" % (k, v) for k, v in params.items()]
        args = shlex.split(' '.join(arg_strs))
        return args

    def test_gets_params_from_cli(self):
        args = self.params_to_args(params=self.params)
        base_command.BaseCommand.run(args=args)
        self.assert_is_subdict(self.mocks['BaseCommand']['handle'].call_args[1],
                               self.params)

    def test_gets_params_from_file(self):
        params_dir = tempfile.mkdtemp(prefix='base_command.params.')
        params_file_path = os.path.join(params_dir, 'params.json')
        with open(params_file_path, 'w') as f: json.dump(self.params, f)
        args = self.params_to_args(params={'params_file': params_file_path})
        base_command.BaseCommand.run(args=args)
        self.assert_is_subdict(self.mocks['BaseCommand']['handle'].call_args[1],
                               self.params)

class CreateFlowRecordCommandTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['odyssey_push_runner'] = patch.multiple(
            odyssey_push_runner, OdysseyPushRunner=DEFAULT)

    def test_handle_calls_runner_create_flow_record(self):
        mock_runner = self.mocks['odyssey_push_runner']['OdysseyPushRunner']\
                .return_value
        mock_flow = {'mock': 'flow'}
        mock_runner.create_flow_record.return_value = mock_flow
        flow_spec_json = json.dumps({'mock': 'flow'})
        command = create_flow_record_command.Command()
        stdout = io.StringIO()
        command.set_streams(stdout=stdout)
        command.handle(flow_spec_json=flow_spec_json)
        self.assertEqual(mock_runner.create_flow_record.call_args,
                         call(flow_record={'spec': flow_spec_json}))
        expected_output = {"status": "SUCCESS",
                           "flow": mock_runner.create_flow_record.return_value}
        self.assertEqual(json.loads(stdout.getvalue().strip()), expected_output)

class RunnerCreateFlowRecordTestCase(BaseTestCase):
    def test_wraps_mc_client_method(self):
        mock_mc_client = Mock()
        self.runner.mc_client = mock_mc_client
        mock_flow_record = Mock()
        result = self.runner.create_flow_record(flow_record=mock_flow_record)
        self.assertEqual(mock_mc_client.create_flow.call_args,
                         call(flow=mock_flow_record))
        self.assertEqual(result, mock_mc_client.create_flow.return_value)

class RunCommandTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['odyssey_push_runner'] = patch.multiple(
            odyssey_push_runner, OdysseyPushRunner=DEFAULT)

    def test_handle_calls_runner_create_flow(self):
        command = run_command.Command()
        command.handle()
        mock_runner = self.mocks['odyssey_push_runner']['OdysseyPushRunner']\
                .return_value
        self.assertEqual(mock_runner.run.call_args, call())

class RunnerSetupTestCase(BaseTestCase):
    def test_falls_back_to_generate_flow_generator_classes(self):
        self.run_subcomponent_generator_fallback_test('flow_generator_classes')

    def run_subcomponent_generator_fallback_test(self, subcomponent_name=None):
        generator_method_name = 'generate_%s' % subcomponent_name
        mock_generator = Mock()
        setattr(self.mock_runner, generator_method_name, mock_generator)
        self.call_real_setup(**{subcomponent_name: Mock()})
        self.assertEqual(mock_generator.call_count, 0)
        self.call_real_setup(**{subcomponent_name: None})
        self.assertEqual(mock_generator.call_count, 1)

    def call_real_setup(self, *args, **kwargs):
        self.mock_runner.call_real_method('setup', *args, **kwargs)

    def test_falls_back_to_generate_task_handler(self):
        self.run_subcomponent_generator_fallback_test('task_handler')

    def test_falls_back_to_generate_flow_engine(self):
        self.run_subcomponent_generator_fallback_test('flow_engine')

    def test_falls_back_to_generate_mc_client(self):
        self.run_subcomponent_generator_fallback_test('mc_client')

    def test_has_job_submission_factory(self):
        self.run_subcomponent_generator_fallback_test('job_submission_factory')

    def test_has_job_runner(self):
        self.run_subcomponent_generator_fallback_test('job_runner')

    def test_has_tick_ctx(self):
        tick_ctx = Mock()
        self.call_real_setup(tick_ctx=tick_ctx)
        self.assertEqual(self.mock_runner.decorate_tick_ctx.call_args,
                         call(tick_ctx=tick_ctx))

    def test_has_flow_runner(self):
        self.run_subcomponent_generator_fallback_test('flow_runner')

class GenerateFlowGeneratorClassesTestCase(BaseTestCase):
    def test_generates_flow_generator_classes(self):
        expected_flow_generator_classes = set()
        self.assertEqual(self.runner.generate_flow_generator_classes(),
                         expected_flow_generator_classes)

class GenerateFlowEngineTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.task_handler = Mock()
        self.runner.flow_generator_classes = [Mock() for i in range(3)]

    def decorate_patchers(self):
        self.patchers['FlowEngine'] = patch.object(odyssey_push_runner,
                                                   'FlowEngine')

    def test_generates_flow_engine(self):
        flow_engine = self.runner.generate_flow_engine()
        self.assertEqual(self.mocks['FlowEngine'].call_args,
                         call(task_handler=self.runner.task_handler))
        self.assertEqual(flow_engine, self.mocks['FlowEngine'].return_value)

    def test_registers_flow_generator_classes_with_engine(self):
        self.runner.generate_flow_engine()
        call_args_list = self.mocks['FlowEngine'].return_value\
                .register_flow_generator_class.call_args_list
        expected_call_args_list = [
            call(flow_generator_class=flow_generator_class)
            for flow_generator_class in self.runner.flow_generator_classes
        ]
        self.assertEqual(call_args_list, expected_call_args_list)

class GenerateMissionControlClientTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['MissionControlClient'] = patch.object(
            odyssey_push_runner, 'MissionControlClient')

    def test_generates_mc_client(self):
        result = self.mock_runner.call_real_method('generate_mc_client')
        self.assertEqual(result,
                         self.mocks['MissionControlClient'].return_value)
        self.assertEqual(self.mocks['MissionControlClient'].call_args,
                         call(base_url=self.mock_runner.mc_server_url,
                              request_client=self.mock_runner.request_client))

class GenerateJobSubmissionFactoryTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['JobSubmissionFactory'] = patch.object(
            odyssey_push_runner, 'JobSubmissionFactory')

    def test_generates_job_submission_factory(self):
        result = self.mock_runner.call_real_method(
            'generate_job_submission_factory')
        self.assertEqual(result,
                         self.mocks['JobSubmissionFactory'].return_value)

class GenerateJobRunnerTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['JobRunner'] = patch.object(odyssey_push_runner,
                                                  'OdysseyPushJobRunner')

    def test_generates_job_runner(self):
        job_runner_kwargs = {'mock': 'kwargs'}
        job_runner = self.mock_runner.call_real_method(
            'generate_job_runner', **{'job_runner_kwargs': job_runner_kwargs})
        self.assertEqual(job_runner, self.mocks['JobRunner'].return_value)
        self.assertEqual(
            self.mocks['JobRunner'].call_args,
            call(
                task_handler=self.mock_runner.task_handler,
                job_client=self.mock_runner.mc_client,
                job_submission_factory=self.mock_runner.job_submission_factory,
                ssh_client=self.mock_runner.ssh_client,
                **job_runner_kwargs
            )
        )

class DecorateTickCtxTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.orig_tick_ctx = {'k1': 'v1', 'k2': 'v2'}
        self.decorated_tick_ctx = self.mock_runner.call_real_method(
            'decorate_tick_ctx', tick_ctx=self.orig_tick_ctx)

    def test_includes_original_tick_ctx(self):
        self.assert_is_subdict(self.orig_tick_ctx, self.decorated_tick_ctx)

    def test_create_job_wraps_mc_client(self):
        args = Mock()
        self.decorated_tick_ctx['create_job'](args)
        self.assertEqual(self.mock_runner.mc_client.create_job.call_args,
                         call(args))

    def test_get_job_wraps_mc_client(self):
        args = Mock()
        self.decorated_tick_ctx['get_job'](args)
        self.assertEqual(
            self.mock_runner.mc_client.fetch_job_by_uuid.call_args,
            call(args))

    def test_create_flow_wraps_mc_client(self):
        args = Mock()
        self.decorated_tick_ctx['create_flow'](args)
        self.assertEqual(self.mock_runner.mc_client.create_flow.call_args,
                         call(args))

    def test_get_flow_wraps_mc_client(self):
        args = Mock()
        self.decorated_tick_ctx['get_flow'](args)
        self.assertEqual(
            self.mock_runner.mc_client.fetch_flow_by_uuid.call_args,
            call(args))

class GenerateFlowRunnerTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['FlowRunner'] = patch.object(odyssey_push_runner,
                                                  'FlowRunner')
    def test_generates_flow_runner(self):
        flow_runner = self.mock_runner.call_real_method('generate_flow_runner')
        self.assertEqual(flow_runner, self.mocks['FlowRunner'].return_value)
        self.assertEqual(self.mocks['FlowRunner'].call_args,
                         call(flow_client=self.mock_runner.mc_client,
                              flow_engine=self.mock_runner.flow_engine,
                              tick_ctx=self.mock_runner.tick_ctx))

class RunTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.tick = Mock()

    def decorate_patchers(self):
        self.patchers['time'] = patch.multiple('time', sleep=DEFAULT)

    def test_calls_tick_at_interval(self):
        tick_intervals = [1, 2]
        tick_counter = 0
        for tick_interval in tick_intervals:
            self.runner.tick_interval = tick_interval
            self.runner.run(ntimes=1)
            tick_counter += 1
            self.assertEqual(self.runner.tick.call_count, tick_counter)
        self.assertEqual(
            self.mocks['time']['sleep'].call_args_list,
            [call(tick_interval) for tick_interval in tick_intervals])

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.flow_runner = Mock()
        self.runner.job_runner = Mock()

    def test_ticks_flow_runner(self):
        self.assertEqual(self.runner.flow_runner.tick.call_count, 0)
        self.runner.tick()
        self.assertEqual(self.runner.flow_runner.tick.call_count, 1)

    def test_ticks_job_runner(self):
        self.assertEqual(self.runner.job_runner.tick.call_count, 0)
        self.runner.tick()
        self.assertEqual(self.runner.job_runner.tick.call_count, 1)

if __name__ == '__main__':
    unittest.main()
