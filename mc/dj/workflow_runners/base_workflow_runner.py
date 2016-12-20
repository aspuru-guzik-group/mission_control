import json
import logging
import time


class BaseWorkflowRunner(object):
    def __init__(self, workflow_client=None, workflow_engine=None,
                 tick_interval=120, max_workflows_per_tick=3, logger=None):
        self.workflow_client = workflow_client
        self.workflow_engine = workflow_engine
        self.tick_interval = tick_interval
        self.max_workflows_per_tick = max_workflows_per_tick
        self.logger = logger or logging

        self._ticking = False
        self.tick_counter = 0

    def start(self):
        self._ticking = True
        self.run()

    def stop(self):
        self._ticking = False

    def run(self, ntimes=None, tick_interval=None):
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep(tick_interval=tick_interval)
        else:
            while self._ticking:
                self._tick_and_sleep(tick_interval=tick_interval)

    def _tick_and_sleep(self, tick_interval=None):
        if tick_interval is None:
            tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        self.tick_counter += 1
        self.logger.debug('tick #%s' % self.tick_counter)
        processed_counter = 0
        for workflow_record in self.fetch_claimable_workflow_records():
            did_process = self.process_claimable_workflow_record(
                workflow_record=workflow_record)
            if did_process: processed_counter += 1
            if processed_counter > self.max_workflows_per_tick: break

    def fetch_claimable_workflow_records(self):
        self.logger.debug('fetch_claimable_workflow_records')
        return self.workflow_client.fetch_claimable_workflows()

    def process_claimable_workflow_record(self, workflow_record=None):
        self.logger.debug('process_claimable_workflow')
        processed = False
        claimed_record = self.claim_workflow_record(workflow_record)
        if claimed_record:
            try:
                updates = self.tick_workflow_record(
                    workflow_record=claimed_record)
                self.update_workflow_record(
                    workflow_record=claimed_record,
                    updates={**updates, 'claimed': False})
                processed = True
            except Exception as exception:
                self.logger.exception(exception)
                self.update_workflow_record(workflow_record=claimed_record,
                                            updates={'status': 'FAILED'})
        return processed

    def claim_workflow_record(self, workflow_record=None):
        self.logger.debug('claim_workflow')
        claimed_records = self.workflow_client.claim_workflows(
            uuids=[workflow_record['uuid']])
        return claimed_records.get(workflow_record['uuid'], False)

    def tick_workflow_record(self, workflow_record=None):
        json_workflow = workflow_record['serialization']
        serialized_workflow = json.loads(json_workflow)
        workflow = self.workflow_engine.deserialize_workflow(
            serialized_workflow=serialized_workflow)
        self.workflow_engine.tick_workflow(workflow=workflow)
        updated_serialization = self.workflow_engine.serialize_workflow(
            workflow=workflow)
        updates = {
            'serialization': json.dumps(updated_serialization),
            'status': updated_serialization['status'],
        }
        return updates

    def update_workflow_record(self, workflow_record=None, updates=None):
        self.workflow_client.update_workflows(updates_by_uuid={
            workflow_record['uuid']: updates})

