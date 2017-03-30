def build_job_submission(*args, job=None, cfg=None, submission_dir=None, **kwargs):
    builder = get_submission_builder(
        job=job, cfg=cfg, submission_dir=submission_dir)
    return builder.build_submission()

def get_submission_builder(*args, **kwargs):
    from .submission_builder import SubmissionBuilder
    return SubmissionBuilder(*args, **kwargs)

def run_job_submission(*args, job=None, cfg=None, submission=None, **kwargs):
    runner = get_submission_runner(job=job, cfg=cfg, submission=submission)
    return runner.run_submission()

def get_submission_runner(*args, **kwargs):
    from .submission_runner import SubmissionRunner
    return SubmissionRunner(*args, **kwargs)



