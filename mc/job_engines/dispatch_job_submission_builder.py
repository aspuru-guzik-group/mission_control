from .base_job_submission_builder import BaseJobSubmissionBuilder


class DispatchJobSubmissionBuilder(BaseJobSubmissionBuilder):
    """
    JobSubmissionBuilder that dispatches to another submission builder based
    on the value of cfg['job_engine']['target_env']
    """
    class UnknownTargetError(Exception): pass

    def _build_job_submission(self):
        target = self.get_target_from_cfg()
        if target == 'bash': return self.build_job_bash_submission()
        raise self.UnknownTargetError("target: {target}".format(target=target))

    def get_target_from_cfg(self):
        return self.cfg.get('job_engine', {}).get('target_env')
