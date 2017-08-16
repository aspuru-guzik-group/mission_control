from collections import defaultdict
import logging
import time
import traceback


class FlowRunner(object):
    """
    A FlowRunner encapsulates logic related to claiming and ticking flows.
    """

    def __init__(self, flow_record_client=None, flow_engine=None,
                 task_ctx=None, tick_interval=120, max_flows_per_tick=3,
                 logger=None):
        """
        Args:
            flow_record_client (mc.clients.flow_record_client): a client for
                doing flow_record operations.
            flow_engine (mc.flows.flow_engine): flow_engine to use for ticking
                flows.  Default: mc.flows.flow_engine.FlowEngine instance.
            task_ctx (dict, optional): extra ctx to pass to tasks.
            tick_interval (int, optional): how often to run a tick, in seconds.
                Default: 120.
            max_flows_per_tick (int, optional)]: maximum number of flows to
                claim per tick. Default: 3.
        """
        self.logger = logger or logging
        self.flow_record_client = flow_record_client
        self.flow_engine = flow_engine or self.get_default_flow_engine()
        self.task_ctx = self.decorate_task_ctx(task_ctx=task_ctx)
        self.tick_interval = tick_interval
        self.max_flows_per_tick = max_flows_per_tick
        self.tick_counter = 0
        self._ticking = False

    def get_default_flow_engine(self):
        from mc.flows.flow_engine import FlowEngine
        return FlowEngine()

    def decorate_task_ctx(self, task_ctx=None):
        task_ctx = task_ctx or {}
        return {**task_ctx, 'flow_engine': self.flow_engine}

    def run(self, ntimes=None, tick_interval=None):
        """Run indefinitely or for several ticks.

        Args:
            ntimes (int, optional): if specified, run this many ticks. If
                empty, run indefinitely. Default: None.
            tick_interval (int, optional): run with this tick interval.
                Default: self.tick_interval.
        """
        self._ticking = True
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep(tick_interval=tick_interval)
        else:
            while self._ticking:
                self._tick_and_sleep(tick_interval=tick_interval)

    def stop(self):
        """Stop ticking."""
        self._ticking = False

    def _tick_and_sleep(self, tick_interval=None):
        if tick_interval is None:
            tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        """Run one claim-and-tick-flows cycle."""
        self.tick_counter += 1
        self.logger.debug('%s, tick #%s' % (self, self.tick_counter))
        claimed_flow_records = self.claim_flow_records()
        tick_stats = {
            'claimed': len(claimed_flow_records),
            **self.tick_flow_records(flow_records=claimed_flow_records)
        }
        return tick_stats

    def claim_flow_records(self):
        return self.flow_record_client.claim_flow_records()

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
        flow = self.flow_record_to_flow(flow_record=flow_record)
        flow.data.setdefault('_flow_record_tick_counter', 0)
        flow.data['_flow_record_tick_counter'] += 1
        self.flow_engine.tick_flow_until_has_no_pending(
            flow=flow, task_ctx=self.task_ctx)
        updated_flow_dict = self.flow_engine.flow_to_flow_dict(flow=flow)
        return updated_flow_dict

    def flow_record_to_flow(self, flow_record=None):
        return self.flow_engine.flow_dict_to_flow(flow_dict=flow_record)

    def patch_and_release_flow_record(self, flow_record=None, patches=None):
        self.flow_record_client.patch_and_release_flow_record(
            flow_record=flow_record, patches=patches)
