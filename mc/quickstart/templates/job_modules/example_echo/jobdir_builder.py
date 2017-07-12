from mc.job_module_utils.jobdir_builder import JobdirBuilder

class ExampleEchoJobdirBuilder(JobdirBuilder):
    def _decorate_job_spec(self, job_spec=None):
        return job_spec

build_jobdir = ExampleEchoJobdirBuilder.build_jobdir
