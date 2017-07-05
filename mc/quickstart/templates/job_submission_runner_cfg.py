def get_cfg_for_job(job=None):
    cfg = {}
    if job['job_spec']['job_type'] == 'job_modules.example_echo':
        cfg = {
            'example_echo': {
                'env_vars': {
                    'ECHO_EXE': 'echo foo'
                }
            }
        }
    return cfg
