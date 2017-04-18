from ...a2g2_common.odyssey import (
    submission_builder as odyssey_submission_builder)
build_job_submission = odyssey_submission_builder.build_job_submission
assert build_job_submission

from . import submission_runner as flow_submission_runner
run_job_submission = flow_submission_runner.run_job_submission
assert run_job_submission
