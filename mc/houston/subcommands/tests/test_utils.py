import unittest
from unittest.mock import call, MagicMock, patch

from .. import _utils


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.cfg = MagicMock()
        self.utils = _utils.HoustonSubcommandUtils(
            get_cfg=MagicMock(return_value=self.cfg)
        )

    def mockify_utils_attrs(self, attrs=None):
        for attr in attrs: setattr(self.utils, attr, MagicMock())

    def mockify_module_attrs(self, attrs=None, module=_utils):
        mod_mocks = {}
        for attr in attrs:
            patcher = patch.object(module, attr)
            self.addCleanup(patcher.stop)
            mod_mocks[attr] = patcher.start()
        return mod_mocks

class McDaoTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mod_mocks = self.mockify_module_attrs(attrs=['_McSqlAlchemyDao'])
        self.result = self.utils.mc_dao

    def test_constructs_dao(self):
        self.assertEqual(self.mod_mocks['_McSqlAlchemyDao'].call_args,
                         call(db_uri=self.cfg['MC_DB_URI']))

    def test_ensures_tables(self):
        self.assertEqual(
            (self.mod_mocks['_McSqlAlchemyDao'].return_value.ensure_tables
             .call_args),
            call()
        )

    def test_returns_dao(self):
        self.assertEqual(self.result,
                         self.mod_mocks['_McSqlAlchemyDao'].return_value)

    def test_memoizes(self):
        self.assertTrue(self.utils.mc_dao is self.result)

class FlowRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_utils_attrs(attrs=['flow_record_client',
                                        'job_record_client'])
        self.mod_mocks = self.mockify_module_attrs(attrs=['FlowRunner'])
        self.result = self.utils.flow_runner

    def test_constructs_and_returns_flow_runner(self):
        self.assertEqual(
            self.mod_mocks['FlowRunner'].call_args,
            call(
                flow_record_client=self.utils.flow_record_client,
                task_ctx={
                    'mc.flow_record_client': self.utils.flow_record_client,
                    'mc.job_record_client': self.utils.job_record_client,
                }
            )
        )
        self.assertEqual(self.result, self.mod_mocks['FlowRunner'].return_value)

    def test_memoizes(self):
        self.assertTrue(self.utils.flow_runner is self.result)

class FlowRecordClientTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_utils_attrs(attrs=['mc_dao'])
        self.mod_mocks = self.mockify_module_attrs(attrs=['FlowRecordClient'])
        self.result = self.utils.flow_record_client

    def test_constructs_and_returns_flow_record_client(self):
        self.assertEqual(
            self.mod_mocks['FlowRecordClient'].call_args,
            call(
                mc_dao=self.utils.mc_dao,
                use_locks=self.utils.cfg.get('USE_LOCKS', True),
                queue_key=self.utils.cfg['FLOW_QUEUE_KEY']
            )
        )
        self.assertEqual(self.result,
                         self.mod_mocks['FlowRecordClient'].return_value)

    def test_memoizes(self):
        self.assertTrue(self.utils.flow_record_client is self.result)

class JobRecordClientTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_utils_attrs(attrs=['mc_dao'])
        self.mod_mocks = self.mockify_module_attrs(attrs=['JobRecordClient'])
        self.result = self.utils.job_record_client

    def test_constructs_and_returns_job_record_client(self):
        self.assertEqual(
            self.mod_mocks['JobRecordClient'].call_args,
            call(
                mc_dao=self.utils.mc_dao,
                use_locks=self.utils.cfg.get('USE_LOCKS', True),
                queue_key=self.utils.cfg['JOB_QUEUE_KEY']
            )
        )
        self.assertEqual(self.result,
                         self.mod_mocks['JobRecordClient'].return_value)

    def test_memoizes(self):
        self.assertTrue(self.utils.job_record_client is self.result)

class JobRunnerTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_utils_attrs(attrs=['job_record_client', 'jobman',
                                        'jobdir_factory'])
        self.mod_mocks = self.mockify_module_attrs(attrs=['JobRunner'])
        self.result = self.utils.job_runner

    def test_constructs_and_returns_job_runner(self):
        self.assertEqual(
            self.mod_mocks['JobRunner'].call_args,
            call(
                artifact_handler=self.utils.cfg['ARTIFACT_HANDLER'],
                job_record_client=self.utils.job_record_client,
                jobman=self.utils.jobman,
                jobdirs_dir=self.utils.cfg.get('JOBDIRS_DIR', None),
                build_jobdir_fn=self.utils.build_jobdir
            )
        )
        self.assertEqual(self.result, self.mod_mocks['JobRunner'].return_value)

    def test_memoizes(self):
        self.assertTrue(self.utils.job_runner is self.result)

class JobManTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mod_mocks = self.mockify_module_attrs(attrs=['JobMan',
                                                          'import_utils'])
        self.result = self.utils.jobman

    def test_loads_jobman_cfg(self):
        self.assertEqual(
            self.mod_mocks['import_utils'].load_module_from_path.call_args,
            call(path=self.utils.cfg['JOBMAN_CFG_PATH'])
        )

    def test_constructs_and_returns_jobman(self):
        self.assertEqual(
            self.mod_mocks['JobMan'].from_cfg.call_args,
            call(cfg=(self.mod_mocks['import_utils'].load_module_from_path
                      .return_value))
        )
        self.assertEqual(self.result,
                         self.mod_mocks['JobMan'].from_cfg.return_value)

    def test_memoizes(self):
        self.assertTrue(self.utils.jobman is self.result)

class BuildJobdirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mod_mocks = self.mockify_module_attrs(attrs=[
            'JobModuleCommandDispatcher'])
        self.args = [MagicMock() for i in range(3)]
        self.kwargs = {'kwarg_%s' % i: MagicMock() for i in range(3)}
        self.result = self.utils.build_jobdir(*self.args, **self.kwargs)

    def test_dispatches_to_job_module_command_dispatcher(self):
        self.assertEqual(self.mod_mocks['JobModuleCommandDispatcher'].call_args,
                         call())
        self.assertEqual(
            (self.mod_mocks['JobModuleCommandDispatcher'].return_value
             .build_jobdir.call_args),
            call(*self.args, **self.kwargs)
        )
        self.assertEqual(
            self.result, 
            (self.mod_mocks['JobModuleCommandDispatcher'].return_value
             .build_jobdir.return_value)
        )
