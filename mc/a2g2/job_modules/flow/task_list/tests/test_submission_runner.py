import collections
import unittest
from unittest.mock import call, MagicMock

from .. import submission_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.job = MagicMock()
        self.cfg = MagicMock()
        self.submission = collections.defaultdict(MagicMock)
        self.submission_runner = submission_runner.SubmissionRunner(
            job=self.job, cfg=self.cfg, submission=self.submission)

class RunSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        for method_name in ['generate_flow', 'generate_flow_ctx', 'run_flow']:
            setattr(self.submission_runner, method_name, MagicMock())
        self.submission_runner.run_submission()

    def test_runs_flow(self):
        self.assertEqual(
            self.submission_runner.run_flow.call_args,
            call(flow=self.submission_runner.generate_flow.return_value,
                 flow_ctx=self.submission_runner.generate_flow_ctx.return_value)
        )

class GenerateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_runner.flow_engine = MagicMock()

    def test_dispatches_to_flow_engine(self):
        flow = self.submission_runner.generate_flow()
        self.assertEqual(
            self.submission_runner.flow_engine.generate_flow.call_args,
            call(flow_spec=\
                 self.submission_runner.job['job_params']['flow_spec']))
        self.assertEqual(
            flow,
            self.submission_runner.flow_engine.generate_flow.return_value)

class GenerateFlowCtxTestCase(BaseTestCase):
    def test_generates_expected_flow_ctx(self):
        flow_ctx = self.submission_runner.generate_flow_ctx()
        expected_flow_ctx = {
            'job': self.job,
            'cfg': self.cfg,
            'submission': self.submission
        }
        self.assertEqual(flow_ctx, expected_flow_ctx)

class RunFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.submission_runner.flow_engine = MagicMock()
        self.flow = MagicMock()
        self.flow_ctx = MagicMock()

    def test_dispatches_to_flow_engine(self):
        self.submission_runner.run_flow(flow=self.flow, flow_ctx=self.flow_ctx)
        self.assertEqual(self.submission_runner.flow_engine.run_flow.call_args,
                         call(flow=self.flow, flow_ctx=self.flow_ctx))

    def test_raises_error_if_flow_fails(self):
        self.fail()
