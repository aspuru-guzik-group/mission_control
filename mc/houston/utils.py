from jobman.jobman import JobMan

from mc.clients.job_record_client import JobRecordClient
from mc.clients.flow_record_client import FlowRecordClient
from mc.flows.flow_engine import FlowEngine
from mc.db.db import Db
from mc.runners.flow_runner import FlowRunner
from mc.runners.jobman_job_runner.job_runner import JobRunner


class HoustonUtils(object):
    def __init__(self, houston=None):
        self.houston = houston

    @property
    def cfg(self): return self.houston.cfg

    @property
    def db(self):
        if not hasattr(self, '_db'):
            self._db = self.generate_db(db_uri=self.cfg['MC_DB_URI'])
        return self._db

    def generate_db(self, db_uri=None, schema=None):
        return Db(db_uri=db_uri, schema=schema)

    @db.setter
    def db(self, value): self._subcommands = value

    def ensure_queues(self):
        self.ensure_queue(queue_cfg=self.cfg['FLOW_QUEUE'])
        self.ensure_queue(queue_cfg=self.cfg['JOB_QUEUE'])

    def ensure_queue(self, queue_cfg=None):
        try:
            self.db.get_item_by_key(item_type='queue', key=queue_cfg['key'])
        except self.db.ItemNotFoundError:
            self.db.create_item(
                item_type='queue',
                item_kwargs={
                    'key': queue_cfg['key'],
                    **queue_cfg.get('queue_kwargs', {})
                }
            )

    @property
    def flow_runner(self):
        if not hasattr(self, '_flow_runner'):
            self._flow_runner = FlowRunner(
                flow_engine=self.flow_engine,
                flow_record_client=self.flow_record_client,
                task_ctx={
                    'mc.flow_record_client': self.flow_record_client,
                    'mc.job_record_client': self.job_record_client,
                }
            )
        return self._flow_runner

    @flow_runner.setter
    def flow_runner(self, new_value): self._flow_runner = new_value

    @property
    def flow_engine(self):
        if not hasattr(self, '_flow_engine'):
            self._flow_engine = FlowEngine()
        return self._flow_engine

    @flow_engine.setter
    def flow_engine(self, new_value): self._flow_engine = new_value

    @property
    def flow_record_client(self):
        if not hasattr(self, '_flow_record_client'):
            self._flow_record_client = self._get_mc_client(record_type='flow')
        return self._flow_record_client

    @flow_record_client.setter
    def flow_record_client(self, new_value):
        self._flow_record_client = new_value

    @property
    def job_record_client(self):
        if not hasattr(self, '_job_record_client'):
            self._job_record_client = self._get_mc_client(record_type='job')
        return self._job_record_client

    def _get_mc_client(self, record_type=None):
        client_cls = None
        if record_type == 'flow':
            client_cls = FlowRecordClient
        elif record_type == 'job':
            client_cls = JobRecordClient
        assert client_cls is not None
        queue_cfg = self.cfg[record_type.upper() + '_QUEUE']
        return client_cls(mc_db=self.db,
                          use_locks=self.cfg.get('USE_LOCKS', True),
                          queue_key=queue_cfg['key'])

    @job_record_client.setter
    def job_record_client(self, new_value): self._job_record_client = new_value

    @property
    def job_runner(self, mc_clients=None):
        if not hasattr(self, '_job_runner'):
            self._job_runner = JobRunner(
                artifact_handler=self.cfg['ARTIFACT_HANDLER'],
                job_record_client=self.job_record_client,
                jobman=self.jobman,
                jobdirs_dir=self.cfg.get('JOBDIRS_DIR', None),
                build_jobdir_fn=self.build_jobdir,
            )
        return self._job_runner

    @job_runner.setter
    def job_runner(self, new_value): self._job_runner = new_value

    @property
    def jobman(self):
        if not hasattr(self, '_jobman'):
            self._jobman = JobMan.from_cfg(cfg=self.cfg['JOBMAN_CFG'])
        return self._jobman

    @jobman.setter
    def jobman(self, new_value): self._jobman = new_value

    def build_jobdir(self, *args, **kwargs):
        try:
            build_jobdir_fn = self.cfg['BUILD_JOBDIR_FN']
        except:
            def build_jobdir_fn(*args, **kwargs):
                return self.houston.run_command('build_job_dir')
        return build_jobdir_fn(*args, **kwargs)

    def has_unfinished_mc_records(self):
        unfinished_records = self.get_unfinished_mc_records()
        for record_type, records in unfinished_records.items():
            if len(records) > 0:
                return True
        return False

    def get_unfinished_mc_records(self):
        return {
            record_type: self._get_unfinished_mc_items(item_type=record_type)
            for record_type in ['flow', 'job']
        }

    def _get_unfinished_mc_items(self, item_type=None):
        return self.db.get_items(item_type=item_type, query={
            'filters': [
                {'field': 'status', 'op': '! IN',
                 'arg': ['FAILED', 'COMPLETED']}
            ]
        })
