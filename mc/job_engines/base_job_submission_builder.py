class BaseJobSubmissionBuilder(object):

    class UnknownTargetError(Exception): pass

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs or {}
        for attr in ['job', 'cfg', 'output_dir']:
            setattr(self, attr, self.kwargs.get(attr))

    def get_target_env_from_cfg(self):
        return self.cfg.get('job_engine', {}).get('target_env')

    def build_job_submission(self): raise NotImplementedError
