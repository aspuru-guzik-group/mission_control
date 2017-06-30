import os

from mc.daos.sqlalchemy_dao import SqlAlchemyDao as _McSqlAlchemyDao
from mc.utils.commands.subcommand_command import SubcommandCommand
from mc.utils import import_utils


class HoustonCommand(SubcommandCommand):
    subcommands = ['create_flow', 'run_until_completed', 'dump_flows']

    class SettingsError(Exception): pass

    @property
    def settings(self):
        if not hasattr(self, '_houston_settings'):
            self._houston_settings = self._load_houston_settings()
        return self._houston_settings

    def _load_houston_settings(self):
        try: return HoustonSettings()
        except Exception as exc: raise self.SettingsError() from exc

    def create_flow(self, args=None, kwargs=None, unparsed_args=None):
        mc_dao = self._get_mc_dao()
        flow_dict = {}
        flow_record = mc_dao.create_item(item_type='Flow',
                                         item_kwargs=flow_dict)
        print("Created flow_record: ", flow_record)

    def _get_mc_dao(self):
        mc_dao = _McSqlAlchemyDao(db_uri=self.settings['MC_DB_URI'])
        mc_dao.ensure_tables()
        return mc_dao

    def run_until_completed(self, args=None, kwargs=None, unparsed_args=None):
        print("run_until_completed")

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
