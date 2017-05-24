import logging
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from .. import flow_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.runner = flow_runner.FlowRunner(
            flow_client=MagicMock(),
            flow_engine=MagicMock(),
            flow_ctx=MagicMock()
        )

class RunTestCase(BaseTestCase):
    @patch.object(flow_runner, 'time')
    def test_run(self, _time):
        tick_intervals = [1, 2]
        self.runner.tick_interval = 10
        for tick_interval in tick_intervals:
            self.runner.tick_interval = tick_interval
            self.runner.run(ntimes=1)
        self.assertEqual(
            _time.sleep.call_args_list,
            [call(tick_interval) for tick_interval in tick_intervals]
        )

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.runner.claim_flow_records = MagicMock()
        self.runner.tick_flow_records = MagicMock()
        self.runner.tick()

    def test_claims_flow_records(self):
        self.assertEqual(self.runner.claim_flow_records.call_args, call())

    def test_ticks_claimed_flow_records(self):
        self.assertEqual(
            self.runner.tick_flow_records.call_args,
            call(flow_records=self.runner.claim_flow_records.return_value)
        )

class ClaimFlowRecordsTestCase(BaseTestCase):
    def test_dispatches_to_flow_client(self):
        self.runner.claim_flow_records()
        self.assertEqual(self.runner.flow_client.claim_flows.call_args, call())

class TickFlowRecordsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_records = [MagicMock() for i in range(3)]
        self.runner.tick_flow_record = MagicMock()
        self.runner.patch_flow_record = MagicMock()
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        super().tearDown()
        logging.disable(logging.NOTSET)

    def _tick_flow_records(self):
        self.runner.tick_flow_records(flow_records=self.flow_records)

    def test_ticks_flow_records(self):
        self._tick_flow_records()
        self.assertEqual(
            self.runner.tick_flow_record.call_args_list,
            [call(flow_record=flow_record) for flow_record in self.flow_records]
        )

    def test_patches_flow_records(self):
        self._tick_flow_records()
        self.assertEqual(
            self.runner.patch_flow_record.call_args_list,
            [
                call(flow_record=flow_record,
                     patches={'claimed': False,
                              **self.runner.tick_flow_record.return_value})
                for flow_record in self.flow_records
            ]
        )

    @patch.object(flow_runner, 'traceback')
    def test_patches_with_failed_status_for_tick_error(self, mock_traceback):
        self.runner.tick_flow_record.side_effect = Exception()
        self._tick_flow_records()
        self.assertEqual(
            self.runner.patch_flow_record.call_args_list,
            [
                call(flow_record=flow_record,
                     patches={'claimed': False,
                              'status': 'FAILED',
                              'error': mock_traceback.format_exc.return_value})
                for flow_record in self.flow_records
            ]
        )

class TickFlowRecordTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_record = MagicMock()
        self.flow = MagicMock()
        self.runner.get_flow_for_flow_record = MagicMock(return_value=self.flow)

    def decorate_patchers(self):
        self.patchers['json'] = patch.multiple('json', dumps=DEFAULT)

    def test_calls_engine_tick_with_get_flow_result(self):
        self.runner.tick_flow_record(self.flow_record)
        expected_flow = self.runner.get_flow_for_flow_record.return_value
        self.assertEqual(
            self.runner.flow_engine.tick_flow.call_args,
            call(flow=expected_flow, flow_ctx=self.runner.flow_ctx))

    def test_includes_status_in_return(self):
        serialization = {'status': 'COMPLETED'}
        self.runner.flow_engine.serialize_flow.return_value = serialization
        result = self.runner.tick_flow_record(self.flow_record)
        self.assertEqual(result['status'], self.flow.status)

    def test_includes_serialization_in_return(self):
        result = self.runner.tick_flow_record(self.flow_record)
        self.assertEqual(result['serialization'],
                         self.runner.flow_engine.serialize_flow.return_value)

class GetFlowForFlowRecordTestCase(BaseTestCase):
    def test_deserializes(self):
        flow_record = {'serialization': MagicMock()}
        flow = self.runner.get_flow_for_flow_record(flow_record=flow_record)
        self.assertEqual(flow,
                         self.runner.flow_engine.deserialize_flow.return_value)
        
class PatchFlowRecordTestCase(BaseTestCase):
    def test_patch_flow_record(self):
        flow_record = MagicMock()
        patches = MagicMock()
        self.runner.patch_flow_record(flow_record=flow_record, patches=patches)
        self.assertEqual(self.runner.flow_client.patch_flow.call_args, 
                         call(key=flow_record['key'], patches=patches))

if __name__ == '__main__':
    unittest.main()
