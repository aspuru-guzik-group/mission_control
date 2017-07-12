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

JOBMAN_CFG_PATH = os.path.join(CFG_DIR, 'jobman_cfg.py')

JOBDIRS_DIR = os.path.join(ROOT_DIR, 'jobdirs')
os.makedirs(JOBDIRS_DIR, exist_ok=True)

JOB_SUBMISSION_RUNNER_EXE = os.path.join(ROOT_DIR, 'run_job_submission.sh')
SUBMISSION_BUILD_TARGET = 'bash'

from mc.artifact_processors.local_path_artifact_processor import (
    LocalPathArtifactProcessor)
ARTIFACT_HANDLER = LocalPathArtifactProcessor()
