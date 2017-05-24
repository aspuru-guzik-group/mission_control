import logging
import time
import traceback

class FlowRunner(object):
    def __init__(self, flow_client=None, flow_engine=None, task_handler=None,
                 flow_ctx=None, tick_interval=120, max_flows_per_tick=3,
                 logger=None):
        self.logger = logger or logging
        self.flow_client = flow_client
        self.flow_engine = flow_engine or self.get_default_flow_engine(
            task_handler=task_handler)
        self.flow_ctx = flow_ctx
        self.tick_interval = tick_interval
        self.max_flows_per_tick = max_flows_per_tick
        self.tick_counter = 0
        self._ticking = False

    def get_default_flow_engine(self, task_handler=None):
        task_handler = task_handler or self.get_default_task_handler()
        from mc.flow_engines.flow_engine import FlowEngine
        return FlowEngine(task_handler=self.task_handler)

    def get_default_task_handler(self):
        from mc.task_handlers.module_dispatch_task_handler import \
                ModuleDispatchTaskHandler
        return ModuleDispatchTaskHandler()

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
        self.tick_flow_records(flow_records=self.claim_flow_records())

    def claim_flow_records(self): return self.flow_client.claim_flows()

    def tick_flow_records(self, flow_records=None):
        for flow_record in flow_records:
            try:
                patches = self.tick_flow_record(flow_record=flow_record)
                self.patch_and_release_flow_record(flow_record=flow_record,
                                                   patches=patches)
            except Exception as exception:
                self.logger.exception(exception)
                failure_patches = {'status': 'FAILED',
                                   'error': traceback.format_exc()}
                self.patch_and_release_flow_record(flow_record=flow_record,
                                                   patches=failure_patches)

    def tick_flow_record(self, flow_record=None):
        self.logger.debug('tick_flow_record')
        flow = self.get_flow_for_flow_record(flow_record=flow_record)
        self.flow_engine.tick_flow(flow=flow, flow_ctx=self.flow_ctx)
        updated_serialization = self.flow_engine.serialize_flow(flow=flow)
        patches = {'serialization': updated_serialization,
                   'status': flow.status}
        return patches

    def get_flow_for_flow_record(self, flow_record=None):
        return self.flow_engine.deserialize_flow(flow_record['serialization'])

    def patch_and_release_flow_record(self, flow_record=None, patches=None):
        self.flow_client.patch_and_release_flow(flow=flow_record,
                                                patches=patches)
