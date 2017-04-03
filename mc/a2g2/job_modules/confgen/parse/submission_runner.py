from ...a2g2_common import base_transform_submission_runner
from . import confgen_parser
from .. import constants

def run_job_submission(*args, **kwargs):
    return base_transform_submission_runner.run_job_submission(
        *args,
        transform_fn=confgen_parser.parse_confgen_dir,
        output_key=constants.PARSE_OUTPUTS_KEY,
        **kwargs)
