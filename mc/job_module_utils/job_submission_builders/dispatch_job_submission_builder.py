import pprint

from .base_job_submission_builder import BaseJobSubmissionBuilder


SUBMISSION_BUILD_TARGET_CFG_KEY = 'SUBMISSION_BUILD_TARGET'

class DispatchJobSubmissionBuilder(BaseJobSubmissionBuilder):
    """
    JobSubmissionBuilder that dispatches to another submission builder based
    on the value of cfg['job_engine']['target_env']
    """

    class UnknownTargetError(Exception):
        def __init__(self, msg=None, target=None, cfg=None, **kwargs):
            if msg: msg += '; '
            else: msg = ''
            msg += (
                "Given target was: '{target}'. Did you specify '{cfg_key}' in"
                " cfg passed to `build_job_submission`? Cfg was:\n{cfg}."
            ).format(
                target=target,
                cfg_key=SUBMISSION_BUILD_TARGET_CFG_KEY,
                cfg=pprint.pformat(cfg)
            )
            super().__init__(msg, **kwargs)

    def _build_job_submission(self):
        target = self.get_target_from_cfg()
        if target == 'bash': return self.build_job_bash_submission()
        raise self.UnknownTargetError(target=target, cfg=self.cfg)

    def get_target_from_cfg(self):
        return self.cfg.get(SUBMISSION_BUILD_TARGET_CFG_KEY)
