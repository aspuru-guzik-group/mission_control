import os
import sys
import mc.job_engines.cli
from mc.job_engines.submission_builders.bash import BashSubmissionBuilder


def build_job_submission(*args, output_dir=None, **kwargs):
    print("build_job_submission")
    submission_meta = BashSubmissionBuilder().build_submission(
        submission_spec={
            'env_vars': [('PYTHONPATH', '$PYTHONPATH:' +':'.join(sys.path))],
        },
        cfg={
            'job_engine': {
                'job_engine_exe': ('python -m {job_engine_cli_module}').format(
                    job_engine_cli_module=mc.job_engines.cli.__name__)
            }
        },
        output_dir=output_dir,
    )
    return submission_meta

def run_job_submission(*args, submission_meta=None, **kwargs):
    for log_name in ['stdout', 'stderr']:
        log_path = os.path.join(
            submission_meta['dir'],
            submission_meta['std_log_file_names'][log_name]
        )
        with open(log_path, 'w') as log:
            log.write("run_job_submission:%s" % log_name)
