import os


def get_cfg_value(self, cfg=None, key=None, default=None):
    if key in os.environ: return os.environ[key]
    return (cfg or {}).get(key, default)
