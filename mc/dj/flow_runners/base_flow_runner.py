import json
import logging
import time


class BaseFlowRunner(object):
    def __init__(self, flow_client=None, flow_engine=None,
                 tick_interval=120, max_flows_per_tick=3, logger=None,
                 tick_ctx=None):
        self.flow_client = flow_client
        self.flow_engine = flow_engine
        self.tick_interval = tick_interval
        self.max_flows_per_tick = max_flows_per_tick
        self.logger = logger or logging
        self.tick_ctx = tick_ctx

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
        for flow_record in self.fetch_tickable_flow_records():
            did_process = self.process_flow_record(flow_record)
            if did_process: processed_counter += 1
            if processed_counter > self.max_flows_per_tick: break

    def fetch_tickable_flow_records(self):
        self.logger.debug('fetch_tickable_flow_records')
        return self.flow_client.fetch_flows(query_params={
            'claimed': False, 'tickable': True})

    def process_flow_record(self, flow_record=None):
        self.logger.debug('process_claimable_flow')
        processed = False
        claimed_record = self.claim_flow_record(flow_record)
        if claimed_record:
            try:
                updates = self.tick_flow_record(
                    flow_record=claimed_record)
                self.update_flow_record(
                    flow_record=claimed_record,
                    updates={**updates, 'claimed': False})
                processed = True
            except Exception as exception:
                self.logger.exception(exception)
                self.update_flow_record(flow_record=claimed_record,
                                            updates={'status': 'FAILED'})
        return processed

    def claim_flow_record(self, flow_record=None):
        self.logger.debug('claim_flow')
        claimed_records = self.flow_client.claim_flows(
            uuids=[flow_record['uuid']])
        return claimed_records.get(flow_record['uuid'], False)

    def tick_flow_record(self, flow_record=None):
        json_flow = flow_record['serialization']
        serialized_flow = json.loads(json_flow)
        flow = self.flow_engine.deserialize_flow(
            serialized_flow=serialized_flow)
        self.flow_engine.tick_flow(flow=flow, ctx=self.tick_ctx)
        updated_serialization = self.flow_engine.serialize_flow(
            flow=flow)
        updates = {'serialization': json.dumps(updated_serialization)}
        serialized_status = updated_serialization['status']
        if serialized_status == 'COMPLETED':
            updates['status'] = serialized_status
        return updates

    def update_flow_record(self, flow_record=None, updates=None):
        self.flow_client.update_flows(updates_by_uuid={
            flow_record['uuid']: updates})

