import logging
import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from ..base_flow_runner import BaseFlowRunner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_client = MagicMock()
        self.flow_engine = MagicMock()
        self.runner = BaseFlowRunner(
            flow_client=self.flow_client,
            flow_engine=self.flow_engine)

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def start_patchers(self, patchers):
        mocks = {key: patcher.start() for key, patcher in patchers.items()}
        return mocks

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

class RunTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'time': patch.multiple('time', sleep=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_run(self):
        tick_intervals = [1, 2]
        self.runner.tick_interval = 10
        for tick_interval in tick_intervals:
            self.runner.tick_interval = tick_interval
            self.runner.run(ntimes=1)
        self.assertEqual(
            self.mocks['time']['sleep'].call_args_list,
            [call(tick_interval) for tick_interval in tick_intervals])

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'runner': patch.multiple(self.runner,
                                     fetch_tickable_flow_records=DEFAULT,
                                     process_flow_record=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_fetches_tickable_flow_records(self):
        self.runner.tick()
        self.assertTrue(self.mocks['runner'].get(
            'fetch_tickable_flow_records').called)

    def test_processes_fetched_records(self):
        self.runner.max_executing_flows = 3
        extra_flows = 2
        flow_records = []
        for i in range(self.runner.max_executing_flows + extra_flows):
            flow_records.append({'uuid': i})
        self.mocks['runner']['fetch_tickable_flow_records'].return_value = \
                flow_records
        self.runner.tick()
        expected_call_args_list = [
            call(flow_record=flow_record)
            for flow_record in flow_records[:-1 * extra_flows]]
        self.assertTrue(
            self.mocks['runner']['process_flow_record'].call_args_list,
            expected_call_args_list)

class FetchTickableFlowRecordsTestCase(BaseTestCase):
    def test_fetch_tickable_flow_records(self):
        self.runner.fetch_tickable_flow_records()
        self.assertEqual(
            self.flow_client.fetch_flows.call_args, 
            call(query_params={'claimed': False, 'tickable': True}))

class ProcessRecordBase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'runner': patch.multiple(
                self.runner, tick_flow_record=DEFAULT,
                update_flow_record=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.flow_record = MagicMock()
        self.claimed = MagicMock()

class ProcessFlowRecordsTestCase(ProcessRecordBase):
    def setUp(self):
        super().setUp()
        self.flow_client.claim_flows.return_value = {
            self.flow_record['uuid']: self.claimed}

    def tearDown(self):
        logging.disable(logging.NOTSET)
        super().tearDown()

    def test_calls_tick_flow_record(self):
        self.runner.process_flow_record(self.flow_record)
        self.assertEqual(self.mocks['runner']['tick_flow_record'].call_args,
                         call(flow_record=self.claimed))

    def test_calls_update(self):
        mock_updates = {'junk': 'updates'}
        self.mocks['runner']['tick_flow_record'].return_value = mock_updates
        self.runner.process_flow_record(self.flow_record)
        self.assertEqual(
            self.mocks['runner']['update_flow_record'].call_args,
            call(flow_record=self.claimed, updates={
                **mock_updates, 'claimed': False}))

    def test_handles_tick_exception(self):
        logging.disable(logging.ERROR) # silently throw exception.
        exception = Exception("some exception")
        self.mocks['runner']['tick_flow_record'].side_effect = exception
        self.runner.process_flow_record(self.flow_record)
        self.assertEqual(
            self.mocks['runner']['update_flow_record'].call_args,
            call(flow_record=self.claimed, updates={'status': 'FAILED'}))

class ProcessUnclaimableFlowRecordsTestCase(ProcessRecordBase):
    def setUp(self):
        super().setUp()
        self.flow_client.claim_flows.return_value = {
            self.flow_record['uuid']: None}

    def test_does_not_call_tick_flow_record(self):
        self.runner.process_flow_record(self.flow_record)
        self.assertEqual(
            self.mocks['runner']['tick_flow_record'].call_count, 0)

class TickFlowRecord(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'json': patch.multiple('json', loads=DEFAULT, dumps=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.flow_record = MagicMock()

    def test_calls_engine_tick_with_deserialization(self):
        self.runner.tick_flow_record(self.flow_record)
        self.assertEqual(self.mocks['json']['loads'].call_args,
                         call(self.flow_record['serialization']))
        self.assertEqual(
            self.flow_engine.deserialize_flow.call_args,
            call(serialized_flow=self.mocks['json']['loads'].return_value))
        expected_deserialized_flow = self.flow_engine\
                .deserialize_flow.return_value
        self.assertEqual(
            self.flow_engine.tick_flow.call_args,
            call(flow=expected_deserialized_flow))

    def test_returns_serialization(self):
        result = self.runner.tick_flow_record(self.flow_record)
        serialization = self.flow_engine.serialize_flow.return_value
        self.assertEqual(
            self.mocks['json']['dumps'].call_args, call(serialization))
        self.assertEqual(
            result['serialization'], self.mocks['json']['dumps'].return_value)
        
    def test_returns_status_if_serialization_status_is_completed(self):
        serialization = {'status': 'COMPLETED'}
        self.flow_engine.serialize_flow.return_value = serialization
        result = self.runner.tick_flow_record(self.flow_record)
        self.assertEqual(result['status'], serialization['status'])

    def test_does_not_status_if_serialization_status_is_not_completed(self):
        result = self.runner.tick_flow_record(self.flow_record)
        serialization = self.flow_engine.serialize_flow.return_value
        serialization['status'] = 'RUNNING'
        self.assertTrue('status' not in result)

class UpdateFlowRecordTestCase(BaseTestCase):
    def test_update_flow_record(self):
        flow_record = MagicMock()
        updates = MagicMock()
        self.runner.update_flow_record(flow_record=flow_record,
                                           updates=updates)
        self.assertEqual(
            self.flow_client.update_flows.call_args, 
            call(updates_by_uuid={flow_record['uuid']: updates}))

if __name__ == '__main__':
    unittest.main()
