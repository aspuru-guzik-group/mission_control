import logging
import time
import traceback


class BaseFlowRunner(object):
    def __init__(self, flow_client=None, flow_engine=None,
                 tick_interval=120, max_flows_per_tick=3, logger=None,
                 flow_ctx=None):
        self.flow_client = flow_client
        self.flow_engine = flow_engine
        self.tick_interval = tick_interval
        self.max_flows_per_tick = max_flows_per_tick
        self.logger = logger or logging
        self.flow_ctx = flow_ctx

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
        if tick_interval is None: tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        self.tick_counter += 1
        self.logger.debug('%s, tick #%s' % (self, self.tick_counter))
        self.tick_flow_records(flow_records=self.claim_flow_records())

    def claim_flow_records(self):
        return self.flow_client.claim_flows()

    def tick_flow_records(self, flow_records=None):
        common_patches = {'claimed': False}
        for flow_record in flow_records:
            try:
                patches = self.tick_flow_record(flow_record=flow_record)
                self.patch_flow_record(
                    flow_record=flow_record,
                    patches={**common_patches, **patches}
                )
            except Exception as exception:
                self.logger.exception(exception)
                self.patch_flow_record(
                    flow_record=flow_record,
                    patches={**common_patches,
                             'status': 'FAILED',
                             'error': traceback.format_exc()}
                )

    def tick_flow_record(self, flow_record=None):
        self.logger.debug('tick_flow_record')
        flow = self.get_flow_for_flow_record(flow_record=flow_record)
        self.flow_engine.tick_flow(flow=flow, flow_ctx=self.flow_ctx)
        updated_serialization = self.flow_engine.serialize_flow(flow=flow)
        patches = {
            'serialization': updated_serialization,
            'status': flow.status
        }
        return patches

    def get_flow_for_flow_record(self, flow_record=None):
        return self.flow_engine.deserialize_flow(flow_record['serialization'])

    def patch_flow_record(self, flow_record=None, patches=None):
        self.flow_client.patch_flow(key=flow_record['uuid'], patches=patches)
