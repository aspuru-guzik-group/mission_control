from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock, patch

from .. import switch_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.control_value = MagicMock()
        self.task_params = defaultdict(MagicMock, **{
            'control_value': self.control_value
        })
        self.task = defaultdict(MagicMock, **{'task_params': self.task_params})
        self.task_ctx = self.generate_task_ctx()
        self.task_handler = self.generate_task_handler()

    def generate_task_ctx(self, task=None):
        task = task or self.task
        return defaultdict(MagicMock, **{'task': self.task})

    def generate_task_handler(self, task_ctx=None):
        task_ctx = task_ctx or self.task_ctx
        return switch_task_handler.SwitchTaskHandler(task_ctx=task_ctx)

class GenerateProxiedTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.get_matching_case = MagicMock()
        self.result = self.task_handler.generate_proxied_task()

    def test_calls_get_matching_case(self):
        self.assertEqual(self.task_handler.get_matching_case.call_args,
                         call())

    def test_returns_task_for_matching_case(self):
        self.assertEqual(
            self.result,
            self.task_handler.get_matching_case.return_value['task']
        )

class GetMatchingCaseTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.cases = [MagicMock() for i in range(3)]
        self.task_handler.task['task_params']['cases'] = self.cases
        self.task_handler.get_control_value = MagicMock()
        self.task_handler.evaluate_case = MagicMock()

    def _get_matching_case(self):
        return self.task_handler.get_matching_case()

    def test_calls_evaluate_case(self):
        self.task_handler.evaluate_case.return_value = False
        expected_control_value = \
                self.task_handler.get_control_value.return_value
        self._get_matching_case()
        self.assertEqual(self.task_handler.evaluate_case.call_args_list,
                         [call(case=case, control_value=expected_control_value)
                          for case in self.cases])

    def test_returns_first_case_that_passes(self):
        self.case_counter = 0
        def mock_evaluate_case(*args, **kwargs):
            self.case_counter += 1
            if self.case_counter == len(self.cases): return True
            return False
        self.task_handler.evaluate_case.side_effect = mock_evaluate_case
        result = self._get_matching_case()
        self.assertEqual(result, self.cases[self.case_counter - 1])

    def test_returns_default_if_no_case_passes(self):
        self.task_handler.evaluate_case.return_value = False
        result = self._get_matching_case()
        self.assertEqual(result, self.task['task_params']['default_case'])

    def test_fails_if_no_matching_case_and_no_default(self):
        self.task['task_params'] = {'control_value': MagicMock()}
        with self.assertRaises(self.task_handler.NoMatchingCaseError):
            self._get_matching_case()

class EvaluateCaseTestCase(BaseTestCase):
    @patch.object(switch_task_handler, 'operator')
    def test_dispatches_to_op_with_control_value_and_arg(self, _operator):
        case = {'condition': defaultdict(MagicMock, **{'op': 'some_op'})}
        control_value = MagicMock()
        result = self.task_handler.evaluate_case(
            case=case, control_value=control_value)
        expected_op = getattr(_operator, case['condition']['op'])
        self.assertEqual(expected_op.call_args,
                         call(control_value, case['condition']['arg']))
        self.assertEqual(result, expected_op.return_value)

class OnProxiedTaskFinishedTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.proxied_task = MagicMock()

    @patch.object(switch_task_handler.BaseProxyingTaskHandler,
                  'on_proxied_task_finished')
    def test_calls_super_method(self, _mock_on_finished):
        self.task_handler.on_proxied_task_finished()
        self.assertEqual(_mock_on_finished.call_args, call())

    def test_propagates_task_data(self):
        self.task_handler.on_proxied_task_finished()
        self.assertEqual(self.task['data'],
                         self.task['proxied_task'].get('data'))

