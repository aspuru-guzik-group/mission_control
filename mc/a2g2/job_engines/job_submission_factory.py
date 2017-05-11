import tempfile

from . import a2g2_job_engine

class JobSubmissionFactory(object):
    def __init__(self, cfg=None):
        self.cfg = cfg

    @classmethod
    def generate_from_spec(cls, spec=None):
        params = spec['params']
        cfg = {
            'job_engine': {
                'engine_class_spec': params.get('engine_class_spec'),
                'entrypoint_preamble': params.get('entrypoint_preamble'),
                'job_engine_exe': params.get('job_engine_exe')
            }
        }
        return cls(cfg=cfg)

    def build_job_submission(self, job=None, submission_dir=None):
        submission_dir = submission_dir or tempfile.mkdtemp(prefix='sf.')
        job_engine = a2g2_job_engine.A2G2JobEngine()
        submission_meta = job_engine.build_job_submission(
            job=job,
            cfg=self.cfg,
            submission_dir=submission_dir
        )
        return submission_meta

    def generate_job_submission_cfg(self):  
        cfg = {'job_engine': {}}
        return cfg
