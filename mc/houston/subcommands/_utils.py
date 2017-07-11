from jobman.jobman import JobMan

from mc.artifact_processors.local_path_artifact_processor import (
    LocalPathArtifactProcessor)
from mc.clients.job_record_client import JobRecordClient
from mc.clients.flow_record_client import FlowRecordClient
from mc.daos.sqlalchemy_dao import SqlAlchemyDao as _McSqlAlchemyDao
from mc.job_module_utils.dispatcher import JobModuleCommandDispatcher
from mc.runners.flow_runner import FlowRunner
from mc.runners.jobman_job_runner.job_runner import JobRunner
from mc.utils import import_utils


class HoustonSubcommandUtils(object):
    def __init__(self, get_cfg=None):
        self.get_cfg = get_cfg

        self._cfg = ...
        self._mc_dao = ...

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

    def get_flow_runner(self, mc_clients=None):
        mc_clients = mc_clients or self.get_mc_clients()
        task_ctx = {
            'mc.flow_record_client': mc_clients['flow'],
            'mc.job_record_client': mc_clients['job'],
        }
        return FlowRunner(flow_record_client=mc_clients['flow'],
                          task_ctx=task_ctx)

    def get_mc_clients(self, mc_dao=None):
        mc_dao = mc_dao or self.mc_dao
        return {
            'flow': self.get_flow_record_client(mc_dao=mc_dao),
            'job': self.get_job_record_client(mc_dao=mc_dao)
        }

    def get_flow_record_client(self, mc_dao=None, queue_key=...):
        mc_dao = mc_dao or self.mc_dao
        if queue_key is ...: queue_key = self.get_mc_queues['flow']['key']
        return FlowRecordClient(mc_dao=mc_dao, use_locks=True,
                                queue_key=queue_key)

    def get_mc_queues(self): raise NotImplementedError

    def _get_job_record_client(self, mc_dao=None, queue_key=...):
        mc_dao = mc_dao or self.mc_dao
        if queue_key is ...: queue_key = self.get_mc_queues['job']['key']
        return JobRecordClient(mc_dao=mc_dao, use_locks=True,
                               queue_key=queue_key)

    def _get_job_runner(self, mc_clients=None):
        return JobRunner(
            artifact_processor=self._get_artifact_processor(),
            job_record_client=mc_clients['job'],
            jobman=self._get_jobman(),
            submissions_dir=self.cfg['SUBMISSIONS_DIR'],
            submission_factory=self._get_submission_factory(),
        )

    def _get_artifact_processor(self): return LocalPathArtifactProcessor()

    def _get_jobman(self):
        jobman_cfg = import_utils.load_module_from_path(
            path=self.cfg['JOBMAN_CFG_PATH'])
        return JobMan.from_cfg(cfg=jobman_cfg)

    def _get_submission_factory(self):
        class MySubmissionFactory(object):
            def __init__(self_, cfg=None):
                self_.cfg = cfg

            def build_job_submission(self_, job=None, output_dir=None):
                return JobModuleCommandDispatcher().build_job_submission(
                    job=job, cfg=self_.cfg, output_dir=output_dir)
        
        return MySubmissionFactory(cfg={
            'JOB_SUBMISSION_RUNNER_EXE': (
                self.cfg['JOB_SUBMISSION_RUNNER_EXE']),
            'SUBMISSION_BUILD_TARGET': 'bash'
        })

    def _has_unfinished_items(self, mc_dao=None):
        for item_type in ['Flow', 'Job']:
            unfinished_items = self._get_unfinished_items(mc_dao=mc_dao,
                                                          item_type=item_type)
            if unfinished_items: return True
        return False

    def _get_unfinished_items(self, mc_dao=None, item_type=None):
        unfinished_filter = {'prop': 'status', 'op': '! IN',
                             'arg': ['COMPLETED', 'FAILED']}
        return mc_dao.get_items(item_type=item_type,
                                query={'filters':  [unfinished_filter]})

    def _get_failed_flows(self, mc_dao=None):
        return self._get_failed_items(mc_dao=mc_dao, item_type='Flow')

    def _get_failed_items(self, mc_dao=None, item_type=None):
        failed_filter = {'prop': 'status', 'op': '=', 'arg': 'FAILED'}
        return mc_dao.get_items(item_type=item_type,
                                query={'filters':  [failed_filter]})
