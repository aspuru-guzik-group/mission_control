import contextlib
import os
import sys


@contextlib.contextmanager
def add_path(path):
    sys.path.insert(0, path)
    yield
    sys.path.remove(path)

with add_path(os.path.dirname(__file__)):
    import houston_cfg

jobman_db_uri = os.path.join(houston_cfg.DATA_DIR, 'jobman.sqlite.db')

from jobman.engines.local_engine import LocalEngine
local_engine_db_uri = os.path.join(houston_cfg.DATA_DIR,
                                   'jobman.local_engine.db.sqlite')
engine = LocalEngine(db_uri=local_engine_db_uri)

job_engine_states_ttl = .1
