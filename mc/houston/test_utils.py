from jobman.engines.local_engine import LocalEngine


def generate_test_cfg():
    cfg = {
        'MC_DB_URI': 'sqlite://',
        'JOBMAN_CFG': {
            'jobman_db_uri': 'sqlite://',
            'engine': LocalEngine(db_uri='sqlite://')
        },
        'FLOW_QUEUE': {
            'key': 'houston_queue_flow',
            'queue_kwargs': {
                'queue_spec': {'item_type': 'flow'}
            }
        },
        'JOB_QUEUE': {
            'key': 'houston_queue_job',
            'queue_kwargs': {
                'queue_spec': {'item_type': 'job'}
            }
        }
    }
    return cfg
