import os

import settings


jobman_db_uri = os.path.join(settings.DATA_DIR, 'jobman.sqlite.db')

from jobman.engines.local_engine import LocalEngine
local_engine_db_uri = os.path.join(settings.DATA_DIR,
                                   'jobman.local_engine.db.sqlite')
engine = LocalEngine(db_uri=local_engine_db_uri)
