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

    def generate_runner(self, *args, **runner_kwargs):
        return odyssey_push_runner.OdysseyPushRunner(**runner_kwargs)

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
        flow_client = Mock()
        self.runner = self.generate_runner(flow_client=flow_client)
        flow = Mock()
        result = self.runner.create_flow(flow=flow)
        self.assertEqual(flow_client.create_flow.call_args, call(flow=flow))
        self.assertEqual(result, flow_client.create_flow.return_value)

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
    def patch_runner(self, **patches):
        patcher = patch.multiple(odyssey_push_runner.OdysseyPushRunner,
                                 **patches)
        self.patchers['runner'] = patcher
        patched = patcher.start()
        return patched

    def run_subcomponent_generator_fallback_test(self, subcomponent_name=None):
        generator_method_name = 'generate_%s' % subcomponent_name
        patched = self.patch_runner(**{generator_method_name: DEFAULT})
        self.generate_runner(**{subcomponent_name: Mock()})
        self.assertEqual(patched[generator_method_name].call_count, 0)
        self.generate_runner(**{subcomponent_name: None})
        self.assertEqual(patched[generator_method_name].call_count, 1)

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

class GenerateFlowEngineTestCase(BaseTestCase):
    def test_generates_flow_engine(self): self.fail()

class GenerateFlowClientTestCase(BaseTestCase):
    def test_generates_flow_client(self): self.fail()

class GenerateJobClientTestCase(BaseTestCase):
    def test_generates_job_client(self): self.fail()

class GenerateJobRunnerTestCase(BaseTestCase):
    def test_generates_job_runner(self): self.fail()

class GenerateFlowRunnerTestCase(BaseTestCase):
    def test_generates_flow_runner(self): self.fail()
    def test_registers_flow_classes(self): self.fail()

class RunTestCase(BaseTestCase):
    def test_calls_tick_at_interval(self):
        self.fail()

class TickTestCase(BaseTestCase):
    def test_ticks_flow_runner(self):
        self.fail()

    def test_ticks_job_runner(self):
        self.fail()

if __name__ == '__main__':
    unittest.main()
