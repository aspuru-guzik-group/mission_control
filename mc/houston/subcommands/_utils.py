from jobman.jobman import JobMan

from mc.clients.job_record_client import JobRecordClient
from mc.clients.flow_record_client import FlowRecordClient
from mc.daos.sqlalchemy_dao import SqlAlchemyDao as _McSqlAlchemyDao
from mc.job_module_utils.job_module_command_dispatcher import (
    JobModuleCommandDispatcher)
from mc.runners.flow_runner import FlowRunner
from mc.runners.jobman_job_runner.job_runner import JobRunner
from mc.utils import import_utils


class HoustonSubcommandUtils(object):
    def __init__(self, load_cfg=None):
        self.load_cfg = load_cfg

        self._cfg = ...
        self._mc_dao = ...
        self._flow_runner = ...
        self._flow_record_client = ...
        self._job_record_client = ...
        self._flow_queue = ...
        self._job_queue = ...
        self._job_runner = ...
        self._jobman = ...

    @property
    def cfg(self):
        if self._cfg is ...: self._cfg = self.load_cfg()
        return self._cfg

    def ensure_db(self): self.mc_dao.ensure_tables()

    def ensure_queues(self):
        self.ensure_queue(queue_cfg=self.cfg['FLOW_QUEUE'])
        self.ensure_queue(queue_cfg=self.cfg['JOB_QUEUE'])

    def ensure_queue(self, queue_cfg=None):
        try: self.mc_dao.get_item_by_key(item_type='Queue',
                                         key=queue_cfg['key'])
        except self.mc_dao.ItemNotFoundError:
            self.mc_dao.create_item(
                item_type='Queue',
                item_kwargs={
                    'key': queue_cfg['key'],
                    **queue_cfg.get('queue_kwargs', {})
                }
            )

    @property
    def mc_dao(self):
        if self._mc_dao is ...: 
            self._mc_dao = _McSqlAlchemyDao(db_uri=self.cfg['MC_DB_URI'])
        return self._mc_dao

    @mc_dao.setter
    def mc_dao(self, new_value): self._mc_dao = new_value

    @property
    def flow_runner(self):
        if self._flow_runner is ...: 
            self._flow_runner = FlowRunner(
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
    def flow_record_client(self):
        if self._flow_record_client is ...:
            self._flow_record_client = self._get_mc_client(record_type='Flow')
        return self._flow_record_client

    @flow_record_client.setter
    def flow_record_client(self, new_value):
        self._flow_record_client = new_value

    @property
    def job_record_client(self):
        if self._job_record_client is ...:
            self._job_record_client = self._get_mc_client(record_type='Job')
        return self._job_record_client

    def _get_mc_client(self, record_type=None):
        client_cls = None
        if record_type == 'Flow': client_cls = FlowRecordClient
        elif record_type == 'Job': client_cls = JobRecordClient
        assert client_cls is not None
        queue_cfg = self.cfg[record_type.upper() + '_QUEUE']
        return client_cls(mc_dao=self.mc_dao,
                          use_locks=self.cfg.get('USE_LOCKS', True),
                          queue_key=queue_cfg['key'])

    @job_record_client.setter
    def job_record_client(self, new_value): self._job_record_client = new_value

    @property
    def job_runner(self, mc_clients=None):
        if self._job_runner is ...:
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
        if self._jobman is ...:
            jobman_cfg = import_utils.load_module_from_path(
                path=self.cfg['JOBMAN_CFG_PATH'])
            self._jobman = JobMan.from_cfg(cfg=jobman_cfg)
        return self._jobman

    @jobman.setter
    def jobman(self, new_value): self._jobman = new_value

    def build_jobdir(self, *args, **kwargs):
        return JobModuleCommandDispatcher().build_jobdir(*args, **kwargs)

    def has_unfinished_mc_records(self):
        unfinished_records = self.get_unfinished_mc_records()
        for record_type, records in unfinished_records.items():
            if len(records) > 0: return True
        return False

    def get_unfinished_mc_records(self):
        return {
            record_type: self._get_unfinished_mc_items(item_type=record_type)
            for record_type in ['Flow', 'Job']
        }

    def _get_unfinished_mc_items(self, item_type=None):
        return self.mc_dao.get_items(item_type=item_type, query={
            'filters': [
                {'prop': 'status', 'op': '! IN',
                 'arg': ['FAILED', 'COMPLETED']}
            ]
        })
