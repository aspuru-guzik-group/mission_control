from mc.job_module_utils.job_submission_builders\
        .dispatch_job_submission_builder import DispatchJobSubmissionBuilder
from mc.submission_builders.bash_submission_builder import BashSubmissionBuilder


class JobSubmissionBuilder(DispatchJobSubmissionBuilder):
    def build_job_bash_submission(self):
        bash_builder = BashSubmissionBuilder()
        submission_meta = bash_builder.build_submission(
            cfg=self.cfg, output_dir=self.output_dir)
        return submission_meta
