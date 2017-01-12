import json
import os
import tempfile
import shlex
import unittest
from unittest.mock import call, DEFAULT, Mock, patch

from django.test import TestCase

from .. import odyssey_push_runner
from ..commands import base as base_command
from ..commands import create_flow as create_flow_command
from ..commands import run as run_command


class BaseTestCase(TestCase):
    def setUp(self):
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = self.start_patchers()
        self.runner = self.generate_runner()

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

class CreateFlowCommandTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['odyssey_push_runner'] = patch.multiple(
            odyssey_push_runner, OdysseyPushRunner=DEFAULT)

    def test_handle_calls_runner_create_flow(self):
        flow_spec = Mock()
        command = create_flow_command.Command()
        command.handle(flow_spec=flow_spec)
        mock_runner = self.mocks['odyssey_push_runner']['OdysseyPushRunner']\
                .return_value
        self.assertEqual(mock_runner.create_flow.call_args,
                         call(flow={'spec': flow_spec}))

class RunnerCreateFlowTestCase(BaseTestCase):
    def test_wraps_flow_client_method(self):
        mock_flow_client = Mock()
        self.runner.flow_client = mock_flow_client
        mock_flow = Mock()
        result = self.runner.create_flow(flow=mock_flow)
        self.assertEqual(mock_flow_client.create_flow.call_args,
                         call(flow=mock_flow))
        self.assertEqual(result, mock_flow_client.create_flow.return_value)

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
        mock_runner = Mock(odyssey_push_runner.OdysseyPushRunner)
        generator_method_name = 'generate_%s' % subcomponent_name
        mock_generator = Mock()
        setattr(mock_runner, generator_method_name, mock_generator)
        self.call_real_setup_with_mock(
            mock=mock_runner, **{subcomponent_name: Mock()})
        self.assertEqual(mock_generator.call_count, 0)
        self.call_real_setup_with_mock(
            mock=mock_runner, **{subcomponent_name: None})
        self.assertEqual(mock_generator.call_count, 1)

    def call_real_setup_with_mock(self, mock=None, *args, **kwargs):
        self.call_real_method_with_mock(mock=mock, method_name='setup',
                                        *args, **kwargs)

    def call_real_method_with_mock(self, mock=None, method_name=None,
                                   *args, **kwargs):
        return getattr(mock.__class__, method_name)(mock, *args, **kwargs)

    def test_falls_back_to_generate_flow_engine(self):
        self.run_subcomponent_generator_fallback_test('flow_engine')

    def test_falls_back_to_generate_flow_client(self):
        self.run_subcomponent_generator_fallback_test('flow_client')

    def test_has_job_client(self):
        self.run_subcomponent_generator_fallback_test('job_client')

    def test_has_job_runner(self):
        self.run_subcomponent_generator_fallback_test('job_runner')

    def test_has_flow_runner(self):
        self.run_subcomponent_generator_fallback_test('flow_runner')

class GenerateFlowGeneratorClassesTestCase(BaseTestCase):
    def test_generates_flow_generator_classes(self):
        runner = self.generate_runner()
        from a2g2.flow_generators import reaxys
        expected_flow_generator_classes = set([reaxys.ReaxysFlowGenerator])
        self.assertEqual(runner.generate_flow_generator_classes(),
                         expected_flow_generator_classes)

class GenerateFlowEngineTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['FlowEngine'] = patch.object(odyssey_push_runner,
                                                   'FlowEngine')

    def test_generates_flow_engine(self):
        self.assertEqual(self.runner.generate_flow_engine(),
                         self.mocks['FlowEngine'].return_value)

    def test_registers_flow_generator_classes_with_engine(self):
        flow_generator_classes = [Mock() for i in range(3)]
        self.runner.generate_flow_engine(
            flow_generator_classes=flow_generator_classes)
        call_args_list = self.mocks['FlowEngine'].return_value\
                .register_flow_generator_class.call_args_list
        expected_call_args_list = [
            call(flow_generator_class=flow_generator_class)
            for flow_generator_class in flow_generator_classes
        ]
        self.assertEqual(call_args_list, expected_call_args_list)

class GenerateFlowClientTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['FlowClient'] = patch.object(odyssey_push_runner,
                                                   'FlowClient')

    def test_generates_flow_client(self):
        self.assertEqual(self.runner.generate_flow_client(),
                         self.mocks['FlowClient'].return_value)

class GenerateJobClientTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['JobClient'] = patch.object(odyssey_push_runner,
                                                  'JobClient')

    def test_generates_job_client(self):
        self.assertEqual(self.runner.generate_job_client(),
                         self.mocks['JobClient'].return_value)

class GenerateJobRunnerTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['JobRunner'] = patch.object(odyssey_push_runner,
                                                  'OdysseyPushJobRunner')

    def test_generates_job_client(self):
        job_client = Mock()
        self.runner.odyssey_suer = Mock()
        self.runner.odyssey_host = Mock()
        job_runner = self.runner.generate_job_runner(job_client=job_client)
        self.assertEqual(job_runner, self.mocks['JobRunner'].return_value)
        self.assertEqual(self.mocks['JobRunner'].call_args,
                         call(job_client=job_client,
                              odyssey_user=self.runner.odyssey_user,
                              odyssey_host=self.runner.odyssey_host))

class GenerateFlowRunnerTestCase(BaseTestCase):
    def decorate_patchers(self):
        self.patchers['FlowRunner'] = patch.object(odyssey_push_runner,
                                                  'FlowRunner')
    def test_generates_flow_runner(self):
        job_client = Mock()
        flow_client = Mock()
        flow_engine = Mock()
        flow_runner = self.runner.generate_flow_runner(flow_client=flow_client,
                                                       job_client=job_client,
                                                       flow_engine=flow_engine)
        self.assertEqual(flow_runner, self.mocks['FlowRunner'].return_value)
        self.assertEqual(self.mocks['FlowRunner'].call_args,
                         call(flow_client=flow_client, job_client=job_client,
                              flow_engine=flow_engine))

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
