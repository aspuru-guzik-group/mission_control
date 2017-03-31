from ...a2g2_common.parse import base_parse_submission_runner
from .confgen_parser import parse_confgen_dir

def run_job_submission(*args, **kwargs):
    return base_parse_submission_runner.run_job_submission(
        *args, parse_dir_fn=parse_confgen_dir, **kwargs)
