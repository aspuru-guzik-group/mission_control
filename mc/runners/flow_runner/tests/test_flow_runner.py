from collections import defaultdict
import logging
import unittest
from unittest.mock import call, MagicMock, patch
from .. import flow_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.runner = flow_runner.FlowRunner(
            flow_record_client=MagicMock(),
            flow_engine=MagicMock(),
            task_ctx=MagicMock()
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
        self.result = self.runner.tick()

    def test_claims_flow_records(self):
        self.assertEqual(self.runner.claim_flow_records.call_args, call())

    def test_ticks_claimed_flow_records(self):
        self.assertEqual(
            self.runner.tick_flow_records.call_args,
            call(flow_records=self.runner.claim_flow_records.return_value)
        )

    def test_returns_tick_stats(self):
        expected_tick_stats = {
            'claimed': len(self.runner.claim_flow_records),
            **self.runner.tick_flow_records.return_value
        }
        self.assertEqual(self.result, expected_tick_stats)


class ClaimFlowRecordsTestCase(BaseTestCase):
    def test_dispatches_to_flow_record_client(self):
        self.runner.claim_flow_records()
        self.assertEqual(
            self.runner.flow_record_client.claim_flow_records.call_args,
            call())


class TickFlowRecordsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_records = [MagicMock() for i in range(3)]
        self.runner.tick_flow_record = MagicMock()
        self.runner.patch_and_release_flow_record = MagicMock()
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        super().tearDown()
        logging.disable(logging.NOTSET)

    def _tick_flow_records(self):
        return self.runner.tick_flow_records(flow_records=self.flow_records)

    def test_ticks_flow_records(self):
        self._tick_flow_records()
        self.assertEqual(
            self.runner.tick_flow_record.call_args_list,
            [call(flow_record=flow_record)
             for flow_record in self.flow_records]
        )

    def test_patches_and_releases_flow_records(self):
        self._tick_flow_records()
        self.assertEqual(
            self.runner.patch_and_release_flow_record.call_args_list,
            [
                call(flow_record=flow_record,
                     patches=self.runner.tick_flow_record.return_value)
                for flow_record in self.flow_records
            ]
        )

    def test_returns_tick_stats(self):
        result = self._tick_flow_records()
        expected_tick_stats = defaultdict(int)
        for i, flow_record in enumerate(self.flow_records):
            expected_status = \
                    self.runner.tick_flow_record.return_value.get('status')
            expected_tick_stats[expected_status] += 1
        self.assertEqual(result, expected_tick_stats)

    @patch.object(flow_runner, 'traceback')
    def test_patches_with_failed_status_for_tick_error(self, mock_traceback):
        self.runner.tick_flow_record.side_effect = Exception()
        self._tick_flow_records()
        self.assertEqual(
            self.runner.patch_and_release_flow_record.call_args_list,
            [
                call(flow_record=flow_record,
                     patches={'status': 'FAILED',
                              'error': mock_traceback.format_exc.return_value})
                for flow_record in self.flow_records
            ]
        )


class TickFlowRecordTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_record = MagicMock()
        self.runner.flow_record_to_flow = MagicMock()
        self.expected_flow = self.runner.flow_record_to_flow.return_value

    def test_calls_engine_tick_with_flow_record_to_flow_result(self):
        self.runner.tick_flow_record(self.flow_record)
        self.assertEqual(self.runner.flow_record_to_flow.call_args,
                         call(flow_record=self.flow_record))
        self.assertEqual(
            self.runner.flow_engine.tick_flow_until_has_no_pending.call_args,
            call(flow=self.expected_flow, task_ctx=self.runner.task_ctx))

    def test_returns_flow_dict(self):
        mock_flow_dict = {'key_%s' % i: MagicMock() for i in range(3)}
        self.runner.flow_engine.flow_to_flow_dict.return_value = mock_flow_dict
        result = self.runner.tick_flow_record(self.flow_record)
        self.assertEqual(result, mock_flow_dict)


class FlowRecordToFlowTestCase(BaseTestCase):
    def test_to_flow(self):
        flow_record = MagicMock()
        flow = self.runner.flow_record_to_flow(flow_record=flow_record)
        self.assertEqual(
            self.runner.flow_engine.flow_dict_to_flow.call_args,
            call(flow_dict=flow_record)
        )
        self.assertEqual(
            flow, self.runner.flow_engine.flow_dict_to_flow.return_value)


class PatchAndReleaseFlowRecordTestCase(BaseTestCase):
    def test_patch_flow_record(self):
        flow_record = MagicMock()
        patches = MagicMock()
        self.runner.patch_and_release_flow_record(flow_record=flow_record,
                                                  patches=patches)
        self.assertEqual(
            (self.runner.flow_record_client.patch_and_release_flow_record
             .call_args),
            call(flow_record=flow_record, patches=patches)
        )


if __name__ == '__main__':
    unittest.main()
