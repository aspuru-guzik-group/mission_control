from mc.job_engines.dispatch_job_submission_builder import (
    DispatchJobSubmissionBuilder as _McDispatchBuilder)
from mc.job_engines.submission_builders.bash import BashSubmissionBuilder


class JobSubmissionBuilder(_McDispatchBuilder):
    def build_job_bash_submission(self):
        bash_builder = BashSubmissionBuilder()
        submission_meta = bash_builder.build_submission(
            cfg=self.cfg, output_dir=self.output_dir)
        return submission_meta
