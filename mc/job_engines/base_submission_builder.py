class BaseSubmissionBuilder(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs or {}
        for attr in ['job', 'cfg', 'output_dir']:
            setattr(self, attr, self.kwargs.get(attr))

    def get_target_env_from_cfg(self):
        return self.cfg['job_engine']['target_env']

    def build_submission(self): raise NotImplementedError
