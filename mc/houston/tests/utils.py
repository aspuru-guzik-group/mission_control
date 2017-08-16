from unittest.mock import MagicMock
from jobman.engines.local_engine import LocalEngine

from .. houston import Houston


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
        },
        'ARTIFACT_HANDLER': MagicMock(return_value={}),
        'BUILD_JOBDIR_FN': MagicMock(return_value={}),
    }
    return cfg


def generate_test_houston(cfg_overrides=None, houston_kwargs=None,
                          ensure_tables=True, ensure_queues=True):
    cfg = {**generate_test_cfg(), **(cfg_overrides or {})}
    houston = Houston(cfg=cfg, **(houston_kwargs or {}))
    if ensure_tables:
        houston.db.ensure_tables()
    if ensure_queues:
        houston.utils.ensure_queues()
    return houston
