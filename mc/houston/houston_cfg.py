import os

from mc.utils import import_utils

class HoustonCfg():
    def __init__(self, cfg_path=...):
        self._cfg = {}
        self._sources = [self._cfg]
        if cfg_path is not ...:
            cfg_file_source = self._generate_source_for_file(path=cfg_path)
            self._sources.append(cfg_file_source)
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

    def __setitem__(self, key, val): self._cfg[key] = val
