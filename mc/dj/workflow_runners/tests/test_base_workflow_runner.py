import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch
from ..base_workflow_runner import BaseWorkflowRunner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.workflow_client = MagicMock()
        self.workflow_engine = MagicMock()
        self.runner = BaseWorkflowRunner(
            workflow_client=self.workflow_client,
            workflow_engine=self.workflow_engine)

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
            'runner': patch.multiple(
                self.runner,
                fetch_claimable_workflow_records=DEFAULT,
                process_claimable_workflow_record=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_fetches_claimable_workflow_recordss(self):
        self.runner.tick()
        self.assertTrue(
            self.mocks['runner']['fetch_claimable_workflow_records'].called)

    def test_processes_claimable_workflow_records(self):
        self.runner.max_executing_workflows = 3
        extra_workflows = 2
        workflow_records = []
        for i in range(self.runner.max_executing_workflows + extra_workflows):
            workflow_records.append({'uuid': i})
        self.mocks['runner']['fetch_claimable_workflow_records'].return_value \
                = workflow_records
        self.runner.tick()
        call_args_list = self.mocks['runner'].get(
            'process_claimable_workflow_record').call_args_list,
        expected_call_args_list = [
            call(workflow_record=workflow_record)
            for workflow_record in workflow_records[:-1 * extra_workflows]]
        self.assertTrue(call_args_list, expected_call_args_list)

class FetchClaimableWorkflowRecordsTestCase(BaseTestCase):
    def test_fetch_claimable_workflow_records(self):
        self.runner.fetch_claimable_workflow_records()
        self.assertEqual(
            self.workflow_client.fetch_claimable_workflows.call_count, 1)

class ProcessClaimableRecordsBase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'runner': patch.multiple(
                self.runner, tick_workflow_record=DEFAULT,
                update_workflow_record=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.workflow_record = MagicMock()
        self.claimed = MagicMock()

class ProcessClaimableWorkflowRecordsTestCase(ProcessClaimableRecordsBase):
    def setUp(self):
        super().setUp()
        self.workflow_client.claim_workflows.return_value = {
            self.workflow_record['uuid']: self.claimed}

    def test_calls_tick_workflow_record(self):
        self.runner.process_claimable_workflow_record(self.workflow_record)
        self.assertEqual(self.mocks['runner']['tick_workflow_record'].call_args,
                         call(workflow_record=self.claimed))

    def test_calls_update(self):
        mock_updates = {'junk': 'updates'}
        self.mocks['runner']['tick_workflow_record'].return_value = mock_updates
        self.runner.process_claimable_workflow_record(self.workflow_record)
        self.assertEqual(
            self.mocks['runner']['update_workflow_record'].call_args,
            call(workflow_record=self.claimed, updates={
                **mock_updates, 'claimed': False}))

    def test_handles_tick_exception(self):
        exception = Exception("some exception")
        self.mocks['runner']['tick_workflow_record'].side_effect = exception
        self.runner.process_claimable_workflow_record(self.workflow_record)
        self.assertEqual(
            self.mocks['runner']['update_workflow_record'].call_args,
            call(workflow_record=self.claimed, updates={'status': 'FAILED'}))

class ProcessUnclaimableWorkflowRecordsTestCase(ProcessClaimableRecordsBase):
    def setUp(self):
        super().setUp()
        self.workflow_client.claim_workflows.return_value = {
            self.workflow_record['uuid']: None}

    def test_does_not_call_tick_workflow_record(self):
        self.runner.process_claimable_workflow_record(self.workflow_record)
        self.assertEqual(
            self.mocks['runner']['tick_workflow_record'].call_count, 0)

class TickWorkflowRecord(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'json': patch.multiple('json', loads=DEFAULT, dumps=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.workflow_record = MagicMock()

    def test_calls_engine_tick_with_deserialization(self):
        self.runner.tick_workflow_record(workflow_record=self.workflow_record)
        self.assertEqual(self.mocks['json']['loads'].call_args,
                         call(self.workflow_record['serialization']))
        self.assertEqual(
            self.workflow_engine.deserialize_workflow.call_args,
            call(serialized_workflow=self.mocks['json']['loads'].return_value))
        expected_deserialized_workflow = self.workflow_engine\
                .deserialize_workflow.return_value
        self.assertEqual(
            self.workflow_engine.tick_workflow.call_args,
            call(workflow=expected_deserialized_workflow))

    def test_returns_updates(self):
        result = self.runner.tick_workflow_record(
            workflow_record=self.workflow_record)
        serialization = self.workflow_engine.serialize_workflow.return_value
        self.assertEqual(
            self.mocks['json']['dumps'].call_args,
            call(serialization))
        self.assertEqual(
            result, {
                'serialization': self.mocks['json']['dumps'].return_value,
                'status': serialization['status']})

class UpdateWorkflowRecordTestCase(BaseTestCase):
    def test_update_workflow_record(self):
        workflow_record = MagicMock()
        updates = MagicMock()
        self.runner.update_workflow_record(workflow_record=workflow_record,
                                           updates=updates)
        self.assertEqual(
            self.workflow_client.update_workflows.call_args, 
            call(updates_by_uuid={workflow_record['uuid']: updates}))

if __name__ == '__main__':
    unittest.main()
