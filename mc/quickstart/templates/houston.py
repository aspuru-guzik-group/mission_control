#!/usr/bin/env python

import logging
import os
import time

from jobman.jobman import JobMan

from mc.artifact_processors.local_path_artifact_processor import (
    LocalPathArtifactProcessor)
from mc.clients.job_record_client import JobRecordClient
from mc.clients.flow_record_client import FlowRecordClient
from mc.daos.sqlalchemy_dao import SqlAlchemyDao as _McSqlAlchemyDao
from mc.flows.flow import Flow
from mc.job_engines.job_engine import JobEngine
from mc.runners.flow_runner import FlowRunner
from mc.runners.jobman_job_runner.job_runner import JobRunner
from mc.utils.commands.subcommand_command import SubcommandCommand
from mc.utils import import_utils


class HoustonCommand(SubcommandCommand):
    subcommands = ['ensure_queues', 'create_flow', 'run_until_completed',
                   'dump_flows']

    class SettingsError(Exception): pass

    def __init__(self, *args, logger=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logger or self._get_default_logger()
    
    def _get_default_logger(self):
        logger = logging.getLogger(__name__)
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)
        return logger

    @property
    def settings(self):
        if not hasattr(self, '_houston_settings'):
            self._houston_settings = self._load_houston_settings()
        return self._houston_settings

    def _load_houston_settings(self):
        try: return HoustonSettings()
        except Exception as exc: raise self.SettingsError() from exc

    def ensure_queues(self, args=None, kwargs=None, unparsed_args=None):
        mc_dao = self._get_mc_dao()
        for item_type in ['Flow', 'Job']:
            queue_key = self.settings.get(item_type.upper() + '_QUEUE_KEY')
            try:
                mc_dao.create_item(
                    item_type='Queue',
                    item_kwargs={
                        'key': queue_key,
                        'queue_spec': {'item_type': item_type},
                    }
                )
                self.logger.info("Created {item_type} queue".format(
                    item_type=item_type))
            except mc_dao.IntegrityError:
                self.logger.info(("Queue with key '{queue_key}' already"
                                  " exists.").format(queue_key=queue_key))

    def create_flow(self, args=None, kwargs=None, unparsed_args=None):
        mc_dao = self._get_mc_dao()
        flow_spec = {
            'tasks': [
                {'task_type': 'print', 'task_params': {'msg': i}}
                for i in range(3)
            ]
        }
        flow_dict = Flow.from_flow_spec(flow_spec=flow_spec).to_flow_dict()
        flow_record = mc_dao.create_item(item_type='Flow',
                                         item_kwargs=flow_dict)
        print("Created flow_record: ", flow_record)

    def _get_mc_dao(self):
        mc_dao = _McSqlAlchemyDao(db_uri=self.settings['MC_DB_URI'])
        mc_dao.ensure_tables()
        return mc_dao

    def run_until_completed(self, args=None, kwargs=None, unparsed_args=None):
        max_ticks = 1e3
        tick_interval = .5
        mc_dao = self._get_mc_dao()
        mc_clients = self._get_mc_clients(mc_dao=mc_dao)
        flow_runner = self._get_flow_runner(mc_clients=mc_clients)
        job_runner = self._get_job_runner(mc_clients=mc_clients)
        tick_counter = 0
        while self._has_unfinished_items(mc_dao=mc_dao):
            tick_counter += 1
            self.logger.info('Tick #%s' % (tick_counter))
            if tick_counter > max_ticks: raise Exception("exceed max ticks")
            flow_runner.tick()
            job_runner.tick()
            time.sleep(tick_interval)
        self.logger.info('Completed.')

    def _get_mc_clients(self, mc_dao=None):
        return {
            'flow': self._get_flow_record_client(mc_dao=mc_dao),
            'job': self._get_job_record_client(mc_dao=mc_dao)
        }

    def _get_flow_record_client(self, mc_dao=None):
        return FlowRecordClient(
            mc_dao=mc_dao, use_locks=True,
            queue_key=self.settings['FLOW_QUEUE_KEY']
        )

    def _get_job_record_client(self, mc_dao=None):
        return JobRecordClient(
            mc_dao=mc_dao, use_locks=True,
            queue_key=self.settings['JOB_QUEUE_KEY']
        )

    def _get_flow_runner(self, mc_clients=None):
        task_ctx = {
            'mc.flow_record_client': mc_clients['flow'],
            'mc.job_record_client': mc_clients['job'],

        }
        return FlowRunner(flow_record_client=mc_clients['flow'],
                          task_ctx=task_ctx)

    def _get_job_runner(self, mc_clients=None):
        return JobRunner(
            artifact_processor=self._get_artifact_processor(),
            job_record_client=mc_clients['job'],
            jobman=self._get_jobman(),
            submissions_dir=self.settings['SUBMISSIONS_DIR'],
            submission_factory=self._get_job_engine()
        )

    def _get_artifact_processor(self): return LocalPathArtifactProcessor()

    def _get_jobman(self):
        jobman_cfg = import_utils.load_module_from_path(
            path=self.settings['JOBMAN_CFG_PATH'])
        return JobMan.from_cfg(cfg=jobman_cfg)

    def _get_job_engine(self): return JobEngine()

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

    def dump_flows(self, args=None, kwargs=None, unparsed_args=None):
        if 'keys_to_exclude' not in kwargs:
            kwargs = {**kwargs, 'keys_to_exclude': {'graph'}}
        self._dump_items(item_type='Flow')

    def dump_jobs(self, args=None, kwargs=None, unparsed_args=None):
        if 'keys_to_exclude' not in kwargs:
            kwargs = {**kwargs, 'keys_to_exclude': {'data'}}
        self._dump_items(item_type='Job', **kwargs)

    def _dump_items(self, item_type=None, keys_to_exclude=None, filters=None):
        print('==== ' + item_type.upper() + ' ====')
        mc_dao = self._get_mc_dao()
        keys_to_exclude = keys_to_exclude or {}
        for item in mc_dao.get_items(item_type=item_type):
            if not all([filter_(item) for filter_ in (filters or [])]): continue
            for key, value in item.items():
                if key not in keys_to_exclude:
                    print("{key}: {value}".format(key=key, value=value))
            print('-' * 10)

    def dump_locks(self, args=None, kwargs=None, unparsed_args=None):
        self._dump_items(item_type='Lock', **kwargs)

class HoustonSettings():
    DEFAULT_SETTINGS_FILE_NAME = 'settings.py'

    def __init__(self, settings_file_path=...):
        self._settings = {}
        self._sources = [self._settings]
        if settings_file_path is ...:
            settings_file_path = os.path.join(os.path.dirname(__file__),
                                              self.DEFAULT_SETTINGS_FILE_NAME)
        if settings_file_path:
            settings_file_source = self._generate_source_for_file(
                path=settings_file_path)
            self._sources.append(settings_file_source)
        self._sources.append(os.environ)

    def _generate_source_for_file(self, path=None):
        module = import_utils.load_module_from_path(path=path)
        return module.__dict__

    def get(self, key, default=...):
        try: return self[key]
        except KeyError as exc:
            if 'default' is ...: raise exc
            else: return default

    def __getitem__(self, key):
        for source in self._sources:
            try: return source[key]
            except KeyError: pass
        raise KeyError(key)

    def __setitem__(self, key, val): self._settings[key] = val

if __name__ == '__main__': HoustonCommand.run()
