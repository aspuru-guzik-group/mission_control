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

    def claim_flow_records(self): return self.flow_client.claim_flow_records()

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
        flow = self.flow_engine.flow_dict_to_flow(flow_dict=flow_record)
        flow.data.setdefault('_flow_record_tick_counter', 0)
        flow.data['_flow_record_tick_counter'] += 1
        self.flow_engine.tick_flow_until_has_no_pending(
            flow=flow, task_ctx=self.task_ctx)
        updated_flow_dict = self.flow_engine.flow_to_flow_dict(flow=flow)
        patches = {**updated_flow_dict,
                   'num_tickable_tasks': len(flow.get_tickable_tasks())}
        return patches

    def patch_and_release_flow_record(self, flow_record=None, patches=None):
        self.flow_client.patch_and_release_flow_record(flow_record=flow_record,
                                                       patches=patches)
