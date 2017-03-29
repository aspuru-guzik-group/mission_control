def build_job_submission(*args, job=None, cfg=None, submission_dir=None, **kwargs):
    builder = get_submission_builder(
        job=job, cfg=cfg, submission_dir=submission_dir)
    return builder.build_submission()

def get_submission_builder(job=None, cfg=None, submission_dir=None):
    from .submission_builder import SubmissionBuilder
    return SubmissionBuilder(job=job, cfg=cfg, submission_dir=submission_dir)

def run_job_submission(*args, submission=None, **kwargs):
    runner = get_submission_runner(submission=submission)
    return runner.run_submission()

def get_submission_runner(submission=None):
    from .submission_runner import SubmissionRunner
    return SubmissionRunner(submission=submission)



