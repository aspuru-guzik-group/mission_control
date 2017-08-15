"""Config"""
import os


class Config(object):
    """
    A utility class for getting configs from different sources.
    """
    def __init__(self, cfg=None):
        self.cfg = cfg or {}
        self.sources = [self.cfg, os.environ]

    def get(self, key, default=...):
        try:
            return self[key]
        except KeyError as exc:
            if 'default' is ...:
                raise exc
            else:
                return default

    def __getitem__(self, key):
        for source in self.sources:
            try:
                return self._get_cfg_item_from_source(source, key)
            except KeyError:
                pass
        raise KeyError(key)

    def _get_cfg_item_from_source(self, source, key):
        try:
            return source[key]
        except TypeError:
            try:
                return getattr(source, key)
            except AttributeError:
                raise KeyError(key)

    def __setitem__(self, key, val): self.cfg[key] = val
