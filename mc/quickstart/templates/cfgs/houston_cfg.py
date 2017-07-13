import os


CFG_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(CFG_DIR)
DATA_DIR = os.path.join(ROOT_DIR, 'data')
os.makedirs(DATA_DIR, exist_ok=True)

MC_DB_URI = 'sqlite:///' + os.path.join(DATA_DIR, 'mc.db.sqlite')

FLOW_QUEUE = {
    'key': 'houston_queue_flow',
    'queue_kwargs': {
        'queue_spec': {'item_type': 'Flow'}
    }
}
JOB_QUEUE = {
    'key': 'houston_queue_job',
    'queue_kwargs': {
        'queue_spec': {'item_type': 'Job'}
    }
}

JOBDIRS_DIR = os.path.join(ROOT_DIR, 'jobdirs')
os.makedirs(JOBDIRS_DIR, exist_ok=True)

from mc.artifact_processors.local_path_artifact_processor import (
    LocalPathArtifactProcessor)
ARTIFACT_HANDLER = LocalPathArtifactProcessor()

JOBMAN_CFG = {
    'jobman_db_uri': os.path.join(DATA_DIR, 'jobman.sqlite.db'),
    'job_engine_states_ttl': .1,
}
from jobman.engines.local_engine import LocalEngine
from mc.job_module_utils import cli as _job_module_utils_cli
JOBMAN_CFG['engine'] = LocalEngine(
    db_uri=os.path.join(DATA_DIR, 'jobman.local_engine.db.sqlite'),
    cfg={
        'ENGINE_PREAMBLE': (
            'export PYTHONPATH={ROOT_DIR}:$PYTHONPATH'.format(ROOT_DIR=ROOT_DIR)
        ),
        'MC_RUN_JOBDIR_CMD': (
            'python {job_module_utils_cli_path} run_jobdir'.format(
                job_module_utils_cli_path=_job_module_utils_cli.__file__)
        ),
        'EXAMPLE_ECHO_ECHO_CMD': 'echo with_prefix'
    }
)
