import json
from pathlib import Path

from . import constants


def get_job_dir_component_path(job_dir=None, component_name=None):
    return Path(job_dir, constants.JOB_DIR_COMPONENT_PATHS[component_name])


def load_job_dir_component(job_dir=None, component_name=None):
    component_path = get_job_dir_component_path(job_dir=job_dir,
                                                component_name=component_name)
    with open(component_path) as f:
        if str(component_path).endswith('.json'):
            result = json.load(f)
        else:
            result = f.read()
    return result


def is_job_dir(path=None):
    return Path(path, constants.JOB_DIR_COMPONENT_PATHS['job_meta']).exists()


def evaluate_completed_job_dir(job_dir=None):
    status = error = None
    has_completed_checkpoint = get_job_dir_component_path(
        job_dir=job_dir, component_name='completed_checkpoint').exists()
    if has_completed_checkpoint:
        status = 'COMPLETED'
    else:
        status = 'FAILED'
        try:
            error = load_job_dir_component(
                job_dir=job_dir, component_name='failed_checkpoint')
        except Exception:
            error = '<unknown failure>'
    return {'status': status, 'error': error}
