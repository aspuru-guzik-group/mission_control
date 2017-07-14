from jobman.engines.local_engine import LocalEngine


def generate_test_cfg():
    cfg = {
        'MC_DB_URI': 'sqlite://',
        'JOBMAN_CFG': {
            'jobman_db_uri': 'sqlite://',
            'engine': LocalEngine(db_uri='sqlite://')
        }
    }
    return cfg
