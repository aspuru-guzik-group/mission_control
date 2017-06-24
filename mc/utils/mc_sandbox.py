import logging
import time

from mc.daos.sa_dao import SaDao
from mc.clients.job_record_client import JobRecordClient
from mc.clients.flow_record_client import FlowRecordClient
from mc.flow_engines.flow_engine import FlowEngine
from mc.runners.flow_runner import FlowRunner

class McSandbox(object):
    def __init__(self, mc_db_uri='sqlite://', logger=None):
        self.logger = logger or logging
        self.mc_dao = self.setup_mc_dao(mc_db_uri=mc_db_uri)
        self.queues = self.setup_queues()
        self.flow_record_client = self.setup_flow_record_client()
        self.job_record_client = self.setup_job_record_client()
        self.flow_engine = self.setup_flow_engine()
        self.task_ctx = self.setup_task_ctx()
        self.flow_runner = self.setup_flow_runner()
        self.job_runner = self.setup_job_runner()

    def setup_mc_dao(self, mc_db_uri=None):
        mc_dao = SaDao(db_uri=mc_db_uri)
        mc_dao.ensure_tables()
        return mc_dao

    def setup_queues(self):
        return {
            item_type: self.mc_dao.create_item(
                item_type='Queue',
                kwargs={'queue_spec': {'item_type': item_type}
            })
            for item_type in ['Flow', 'Job']
        }

    def setup_task_ctx(self):
        task_ctx = {
            'mc.job_record_client': self.job_record_client,
            'mc.flow_record_client': self.flow_record_client,
        }
        return task_ctx

    def setup_flow_record_client(self):
        return FlowRecordClient(mc_dao=self.mc_dao,
                                queue_key=self.queues['Flow']['key'])

    def setup_job_record_client(self):
        return JobRecordClient(mc_dao=self.mc_dao,
                               queue_key=self.queues['Job']['key'])

    def setup_flow_engine(self):
        return FlowEngine()

    def setup_flow_runner(self):
        return FlowRunner(flow_engine=self.flow_engine,
                          flow_record_client=self.flow_record_client,
                          task_ctx=self.task_ctx)

    def setup_job_runner(self):
        pass

    def has_incomplete_items(self):
        return self.has_incomplete_flows() or self.has_incomplete_jobs()

    def has_incomplete_flows(self):
        return len(self.get_incomplete_items(item_type='Flow')) > 0

    def get_incomplete_items(self, item_type=None):
        incomplete_filter = {'field': 'status', 'operator': '! IN',
                             'value': ['COMPLETED', 'FAILED']}
        return self.mc_dao.get_items(item_type=item_type,
                                     query={'filters':  [incomplete_filter]})
    def has_incomplete_jobs(self):
        return len(self.get_incomplete_items(item_type='Job')) > 0

    def run_until_completed(self, max_ticks=10, tick_interval=.1,
                            log_ticks=False):
        tick_counter = 0
        while (self.has_incomplete_flows() or self.has_incomplete_jobs()):
            tick_counter += 1
            log_msg = 't{tick_counter}:'.format(tick_counter=tick_counter)
            flow_tick_stats = self.flow_runner.tick()
            if flow_tick_stats['claimed'] > 0: log_msg += 'F'
            job_tick_stats = self.job_runner.tick()
            if job_tick_stats['claimed'] > 0: log_msg += 'J'
            log_msg += ' | '
            if log_ticks: self.logger.warn(log_msg)
            if tick_counter > max_ticks: raise Exception("Exceed max_ticks")
            time.sleep(tick_interval)
