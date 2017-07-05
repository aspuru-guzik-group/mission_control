import os


CFG_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(CFG_DIR)
DATA_DIR = os.path.join(ROOT_DIR, 'data')
os.makedirs(DATA_DIR, exist_ok=True)

MC_DB_URI = 'sqlite:///' + os.path.join(DATA_DIR, 'mc.db.sqlite')

FLOW_QUEUE_KEY = 'houston_queue_flow'
JOB_QUEUE_KEY = 'houston_queue_job'

JOBMAN_CFG_PATH = os.path.join(CFG_DIR, 'jobman_cfg.py')

SUBMISSIONS_DIR = os.path.join(ROOT_DIR, 'submissions')
os.makedirs(SUBMISSIONS_DIR, exist_ok=True)

JOB_SUBMISSION_RUNNER_EXE = os.path.join(ROOT_DIR, 'run_job_submission.sh')
SUBMISSION_BUILD_TARGET = 'bash'
