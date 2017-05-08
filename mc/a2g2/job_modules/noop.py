import time

from .a2g2_common.odyssey import (
    submission_builder as odyssey_submission_builder)
build_job_submission = odyssey_submission_builder.build_job_submission
assert build_job_submission

def run_job_submission(*args, **kwargs): print("NOOP", time.time())
