import unittest
from unittest.mock import call, MagicMock

from ..a2g2_task_handler import A2G2TaskHandler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.task_handler = A2G2TaskHandler()

@unittest.skip
class TickTaskTestCase(BaseTestCase):
    def test_ticks_interpolate_task(self):
        self.fail()

@unittest.skip
class GetHandlerTestCase(BaseTestCase):
    def test_gets_expected_handler(self):
        self.fail()

class InterpolatedTaskParamValueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_context = MagicMock()
        self.task = MagicMock()
        self.task_handler.get_ctx_value = MagicMock()

    def  _interpolate(self, value=None):
        return self.task_handler.interpolate_task_param_value(
            value=value,
            task=self.task,
            task_context=self.task_context)

    def test_does_not_interpolate_plain_values(self):
        orig = 'plain'
        interpolated = self._interpolate(value=orig)
        self.assertEqual(interpolated, orig)

    def test_interpolates_ctx_spec(self):
        ctx_spec = 'some.ctx.spec'
        orig = '_ctx:' + ctx_spec
        interpolated = self._interpolate(value=orig)
        self.assertEqual(self.task_handler.get_ctx_value.call_args,
                         call(ctx=self.task_context, dot_spec=ctx_spec))
        self.assertEqual(interpolated,
                         self.task_handler.get_ctx_value.return_value)

@unittest.skip
class GetCtxValueTestCase(BaseTestCase):
    def test_dispatches_to_dot_spec_loader(self):
        self.fail()
