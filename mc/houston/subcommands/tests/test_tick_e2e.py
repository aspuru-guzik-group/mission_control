import unittest

from mc.houston.tests import utils as _houston_test_utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.houston = _houston_test_utils.generate_test_houston()
        self.common_command_kwargs = {'interval': .01, 'max_ticks': int(1e2)}

    def _generate_flow_spec(self, flow_spec=None):
        flow_spec = flow_spec or {
            'key': 'my_flow',
            'tasks': [
                {
                    'key': 'task_%s' % i,
                    'task_type': 'noop',
                    'task_params': {'msg': 'msg_%s' % i}
                }
                for i in range(3)
            ]
        }
        return flow_spec

    def _create_flow_record(self, flow_spec=None):
        flow_dict = self.houston.utils.flow_engine.flow_spec_to_flow_dict(
            flow_spec=flow_spec)
        flow_record = self.houston.utils.flow_record_client.create_flow_record(
            flow_kwargs=flow_dict)
        return flow_record

    def _create_job_record(self, job_kwargs=None):
        job_kwargs = job_kwargs or {}
        job_record = self.houston.utils.job_record_client.create_job_record(
            job_kwargs=job_kwargs)
        return job_record

    def _run_tick_command(self, **kwargs):
        return self.houston.run_command(
            'tick', **{**self.common_command_kwargs, **kwargs})


class TickFlowSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_spec = self._generate_flow_spec()

    def test_runs_single_tick(self):
        output = self._run_tick_command(tickee='flow',
                                        flow_spec=self.flow_spec)
        expected_flow = self.houston.utils.flow_engine.flow_spec_to_flow(
            flow_spec=self.flow_spec)
        self.houston.utils.flow_engine.tick_flow(flow=expected_flow)
        expected_flow_dict = self.houston.utils.flow_engine.flow_to_flow_dict(
            flow=expected_flow)
        self.assertEqual(output, expected_flow_dict)


class TickFlowRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_spec = self._generate_flow_spec()
        self.flow_record = self._create_flow_record(flow_spec=self.flow_spec)

    def test_runs_for_nticks(self):
        self._run_tick_command(tickee='flow_runner', nticks=3)

    def test_runs_until_finished(self):
        self._run_tick_command(tickee='flow_runner', until_finished=True)


class TickJobRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_runs_for_nticks(self):
        self.job_record = self._create_job_record()
        self._run_tick_command(tickee='job_runner', nticks=3)

    def test_runs_until_finished(self):
        self.skipTest('RETURN TO THIS LATER!')
        self.job_record = self._create_job_record({'status': 'COMPLETED'})
        self._run_tick_command(tickee='job_runner', until_finished=True)


class TickJobManTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_runs_for_nticks(self):
        self._run_tick_command(tickee='jobman', nticks=3)


class TickAllTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_runs_for_nticks(self):
        self._run_tick_command(tickee='all', nticks=3)

    def test_runs_unfil_finished(self):
        self._run_tick_command(tickee='all', until_finished=True)
