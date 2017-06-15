from collections import defaultdict
import logging
import time
import traceback

class FlowRunner(object):
    def __init__(self, flow_client=None, flow_engine=None, task_ctx=None,
                 tick_interval=120, max_flows_per_tick=3, logger=None):
        self.logger = logger or logging
        self.flow_client = flow_client
        self.flow_engine = flow_engine or self.get_default_flow_engine()
        self.task_ctx = task_ctx or {}
        self.fill_in_task_ctx()
        self.tick_interval = tick_interval
        self.max_flows_per_tick = max_flows_per_tick
        self.tick_counter = 0
        self._ticking = False


    def get_default_flow_engine(self, task_handler=None):
        from mc.flow_engines.flow_engine import FlowEngine
        return FlowEngine()

    def fill_in_task_ctx(self):
        self.task_ctx.setdefault('flow_engine', self.flow_engine)

    def run(self, ntimes=None, tick_interval=None):
        self._ticking = True
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep(tick_interval=tick_interval)
        else:
            while self._ticking:
                self._tick_and_sleep(tick_interval=tick_interval)

    def stop(self): self._ticking = False

    def _tick_and_sleep(self, tick_interval=None):
        if tick_interval is None: tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        self.tick_counter += 1
        self.logger.debug('%s, tick #%s' % (self, self.tick_counter))
        claimed_flow_records = self.claim_flow_records()
        tick_stats = {
            'claimed': len(claimed_flow_records),
            **self.tick_flow_records(flow_records=claimed_flow_records)
        }
        return tick_stats

    def claim_flow_records(self): return self.flow_client.claim_flows()

    def tick_flow_records(self, flow_records=None):
        tick_stats = defaultdict(int)
        for flow_record in flow_records:
            try:
                patches = self.tick_flow_record(flow_record=flow_record)
                status = patches.get('status')
            except Exception as exception:
                self.logger.exception(exception)
                status = 'FAILED'
                patches = {'status': status, 'error': traceback.format_exc()}
            self.patch_and_release_flow_record(flow_record=flow_record,
                                               patches=patches)
            tick_stats[status] += 1
        return tick_stats

    def tick_flow_record(self, flow_record=None):
        self.logger.debug('tick_flow_record')
        flow = self.get_flow_for_flow_record(flow_record=flow_record)
        flow.data.setdefault('_tick_counter', 0)
        flow.data['_tick_counter'] += 1
        self.flow_engine.tick_flow(flow=flow, task_ctx=self.task_ctx)
        updated_serialization = self.flow_engine.serialize_flow(flow=flow)
        patches = {'serialization': updated_serialization,
                   'status': flow.status}
        return patches

    def get_flow_for_flow_record(self, flow_record=None):
        return self.flow_engine.deserialize_flow(flow_record['serialization'])

    def patch_and_release_flow_record(self, flow_record=None, patches=None):
        self.flow_client.patch_and_release_flow(flow=flow_record,
                                                patches=patches)


    @staticmethod
    def generate_flow_client_from_mc_dao(mc_dao=None, queue_key=None):
        return McDaoFlowClient(mc_dao=mc_dao, queue_key=queue_key)

class McDaoFlowClient(object):
    def __init__(self, mc_dao=None, queue_key=None):
        self.mc_dao = mc_dao
        self.queue_key = queue_key

    def claim_flows(self):
        return self.mc_dao.claim_queue_items(queue_key=self.queue_key)['items']

    def patch_and_release_flow(self, flow=None, patches=None):
        return self.mc_dao.patch_item(item_type='Flow', key=flow['key'],
                                      patches={'claimed': False, **patches})
