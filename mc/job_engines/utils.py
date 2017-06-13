import json
import os

from . import constants as _job_engine_constants

def read_submission_meta(submission_dir=None, submission_meta_name=None):
    submission_meta_name = submission_meta_name or \
            _job_engine_constants.SUBMISSION_META_NAME
    submission_meta_path = os.path.join(submission_dir, submission_meta_name)
    with open(submission_meta_path) as f: return json.load(f)
