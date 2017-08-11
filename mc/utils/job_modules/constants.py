_prefix = 'A2G2__'

JOB_DIR_COMPONENT_PATHS = {
    'job_key': _prefix + 'JOB_KEY',
    'job_meta': _prefix + 'JOB_META.json',
    'job_spec': 'JOBMAN__JOB_SPEC.json',
    'entrypoint': 'entrypoint.sh',
    'work_dir': 'work_dir',
    'executed_checkpoint': _prefix + 'EXECUTED',
    'completed_checkpoint': _prefix + 'COMPLETED',
    'failure_checkpoint': _prefix + 'FAILED',
    'failed_checkpoint': _prefix + 'FAILED',
    'stdout_log': _prefix + 'STDOUT.log',
    'stderr_log': _prefix + 'STDERR.log',
}
