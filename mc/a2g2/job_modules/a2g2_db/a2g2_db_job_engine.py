from ..utils.job_engine import JobEngine
from . import a2g2_db_uploader


def upload_bulk_files(*args, job=None, cfg=None, **kwargs):
    uploader = a2g2_db_uploader.generate_uploader()
    uploader.upload_bulk_files(**{
        'bulk_file_dir': job['job_spec']['kwargs']['bulk_file_dir'],
        'cfg': cfg,
    })

command_handlers = {
    method.__name__: method
    for method in [upload_bulk_files]
}

def generate_job_engine():
    return JobEngine(command_handlers=command_handlers)
