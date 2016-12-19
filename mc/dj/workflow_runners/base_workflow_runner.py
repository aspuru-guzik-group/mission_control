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
        logging.debug('tick #%s' % self.tick_counter)
        processed_counter = 0
        for workflow_record in self.fetch_claimable_workflow_records():
            did_process = self.process_claimable_workflow_record(
                workflow_record=workflow_record)
            if did_process: processed_counter += 1
            if processed_counter > self.max_workflows_per_tick: break

    def fetch_claimable_workflow_records(self):
        logging.debug('fetch_claimable_workflows')
        return self.workflow_client.fetch_claimable_workflow_records()

    def process_claimable_workflow_record(self, workflow_record=None):
        logging.debug('process_claimable_workflow')
        processed = False
        claimed_record = self.claim_workflow_record(workflow_record)
        if claimed_record:
            try:
                updated_serialization = self.tick_workflow_record(
                    workflow_record=claimed_record)
                self.update_workflow_record(
                    workflow_record=claimed_record,
                    updates={'serialization': updated_serialization})
                processed = True
            except Exception as exception:
                logging.exception(exception)
                self.update_workflow_record(workflow_record=claimed_record,
                                            updates={'status': 'FAILED'})
        return processed

    def claim_workflow_record(self, workflow_record=None):
        logging.debug('claim_workflow')
        claimed_records = self.workflow_client.claim_workflow_records(
            uuids=[workflow_record['uuid']])
        return claimed_records.get(workflow_record['uuid'], False)

    def tick_workflow_record(self, workflow_record=None):
        workflow = self.workflow_engine.deserialize_workflow(
            serialization=workflow_record['serialization'])
        self.workflow_engine.tick_workflow(workflow=workflow)
        updated_serialization = self.workflow_engine.serialize_workflow(
            workflow=workflow)
        return updated_serialization

    def update_workflow_record(self, workflow_record=None, updates=None):
        self.workflow_client.update_workflow_records(updates_by_uuid={
            workflow_record['uuid']: updates})

