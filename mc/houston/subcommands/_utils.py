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
    def __init__(self, get_cfg=None):
        self.get_cfg = get_cfg

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
        if self._cfg is ...: self._cfg = self.get_cfg()
        return self._cfg

    @property
    def mc_dao(self):
        if self._mc_dao is ...: 
            self._mc_dao = _McSqlAlchemyDao(db_uri=self.cfg['MC_DB_URI'])
            self._mc_dao.ensure_tables()
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
            self._flow_record_client = FlowRecordClient(
                mc_dao=self.mc_dao,
                use_locks=self.cfg.get('USE_LOCKS', True),
                queue_key=self.cfg['FLOW_QUEUE_KEY']
            )
        return self._flow_record_client

    @flow_record_client.setter
    def flow_record_client(self, new_value):
        self._flow_record_client = new_value

    @property
    def job_record_client(self):
        if self._job_record_client is ...:
            self._job_record_client = JobRecordClient(
                mc_dao=self.mc_dao,
                use_locks=self.cfg.get('USE_LOCKS', True),
                queue_key=self.cfg['JOB_QUEUE_KEY']
            )
        return self._job_record_client

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
