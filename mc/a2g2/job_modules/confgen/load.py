from ..a2g2_common.odyssey.submission_builder import build_job_submission
assert build_job_submission
from ..a2g2_common import base_transform_submission_runner
from ..a2g2_db import a2g2_db_uploader
from . import constants

def run_job_submission(*args, **kwargs):
    return base_transform_submission_runner.run_job_submission(
        *args,
        transform_fn=a2g2_db_uploader.load_from_input_dir,
        output_key=constants.LOAD_OUTPUTS_KEY,
        **kwargs)
