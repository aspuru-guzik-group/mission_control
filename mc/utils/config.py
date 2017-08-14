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
        try: return self[key]
        except KeyError as exc:
            if 'default' is ...: raise exc
            else: return default

    def __getitem__(self, key):
        for source in self.sources:
            try: return source[key]
            except KeyError: pass
        raise KeyError(key)

    def __setitem__(self, key, val): self.cfg[key] = val
