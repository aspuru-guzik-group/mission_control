from collections import defaultdict
import logging
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from .. import flow_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.runner = flow_runner.FlowRunner(
            flow_client=MagicMock(),
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
    def test_dispatches_to_flow_client(self):
        self.runner.claim_flow_records()
        self.assertEqual(self.runner.flow_client.claim_flows.call_args, call())

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
            [call(flow_record=flow_record) for flow_record in self.flow_records]
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
        self.flow = MagicMock()
        self.runner.get_flow_for_flow_record = MagicMock(return_value=self.flow)

    def decorate_patchers(self):
        self.patchers['json'] = patch.multiple('json', dumps=DEFAULT)

    def test_calls_engine_tick_with_get_flow_result(self):
        self.runner.tick_flow_record(self.flow_record)
        expected_flow = self.runner.get_flow_for_flow_record.return_value
        self.assertEqual(
            self.runner.flow_engine.tick_flow.call_args,
            call(flow=expected_flow, task_ctx=self.runner.task_ctx))

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
        
class PatchAndReleaseFlowRecordTestCase(BaseTestCase):
    def test_patch_flow_record(self):
        flow_record = MagicMock()
        patches = MagicMock()
        self.runner.patch_and_release_flow_record(flow_record=flow_record,
                                                  patches=patches)
        self.assertEqual(
            self.runner.flow_client.patch_and_release_flow.call_args, 
            call(flow=flow_record, patches=patches)
        )

class GenerateFlowClientFromMcDaoTestCase(BaseTestCase):
    @patch.object(flow_runner, 'McDaoFlowClient')
    def test_returns_flow_client_instance(self, _McDaoFlowClient):
        mc_dao = MagicMock()
        queue_key = MagicMock()
        result = flow_runner.FlowRunner.generate_flow_client_from_mc_dao(
            mc_dao=mc_dao, queue_key=queue_key)
        self.assertEqual(_McDaoFlowClient.call_args,
                         call(mc_dao=mc_dao, queue_key=queue_key))
        self.assertEqual(result, _McDaoFlowClient.return_value)

class McDaoFlowClientTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_client = flow_runner.McDaoFlowClient(
            mc_dao=MagicMock(), queue_key=MagicMock())

    def test_client_claim_flows(self):
        result = self.flow_client.claim_flows()
        self.assertEqual(self.flow_client.mc_dao.claim_queue_items.call_args,
                         call(queue_key=self.flow_client.queue_key))
        self.assertEqual(
            result,
            self.flow_client.mc_dao.claim_queue_items.return_value['items'])

    def test_client_patch_and_release_flow(self):
        flow = MagicMock()
        patches = MagicMock()
        result = self.flow_client.patch_and_release_flow(
            flow=flow, patches=patches)
        self.assertEqual(
            self.flow_client.mc_dao.patch_item.call_args,
            call(item_type='Flow', key=flow['key'],
                 patches={'claimed': False, **patches})
        )
        self.assertEqual(result,
                         self.flow_client.mc_dao.patch_item.return_value)

if __name__ == '__main__':
    unittest.main()
