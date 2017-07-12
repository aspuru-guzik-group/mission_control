import contextlib
import os
import sys
import textwrap

from mc.job_module_utils import cli


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
engine = LocalEngine(
    db_uri=local_engine_db_uri,
    cfg={
        'ENGINE_PREAMBLE': textwrap.dedent(
            """
            export PYTHONPATH={project_root_dir}:$PYTHONPATH
            """
        ).lstrip().format(project_root_dir=houston_cfg.ROOT_DIR),
        'MC_RUN_JOBDIR_CMD': (
            'python {job_module_utils_cli_path} run_jobdir'
        ).format(job_module_utils_cli_path=cli.__file__)
    }
)

job_engine_states_ttl = .1
