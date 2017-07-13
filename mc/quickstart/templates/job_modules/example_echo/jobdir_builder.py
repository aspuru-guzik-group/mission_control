from mc.job_module_utils.jobdir_builder import JobdirBuilder

from . import cfg


class ExampleEchoJobdirBuilder(JobdirBuilder):
    def _decorate_job_spec(self, job_spec=None):
        job_spec.setdefault('cfg_specs', {}).update(cfg.cfg_specs)
        return job_spec

build_jobdir = ExampleEchoJobdirBuilder.build_jobdir
