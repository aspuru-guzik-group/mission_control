class BaseJobSubmissionBuilder(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs or {}
        for attr in ['job', 'cfg']:
            setattr(self, attr, self.kwargs.get(attr, {}))
        self.output_dir = self.kwargs.get('output_dir')

    @classmethod
    def build_job_submission(cls, *args, job=None, cfg=None, output_dir=None,
                             **kwargs):
        builder = cls(*args, job=job, cfg=cfg, output_dir=output_dir, **kwargs)
        return builder._build_job_submission()

    def _build_job_submission(self): raise NotImplementedError
    
