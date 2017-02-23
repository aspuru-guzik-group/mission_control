from mc.job_runners.base_job_runner import BaseJobRunner
from mc.job_client.job_client import MissionControlJobClient


def generate_stub_job_runner(base_url=None, request_client=None):
    class StubJobDirFactory(object):
        def build_dir_for_job(self, job=None):
            return {}

    class StubExecutionClient(object):
        def start_execution(self, job=None):
            return {}

        def get_execution_state(self, job=None):
            return {}

    runner = BaseJobRunner(
        execution_client=StubExecutionClient(),
        job_client=MissionControlJobClient(
            base_url=base_url,
            request_client=request_client
        ),
        job_dir_factory=StubJobDirFactory()
    )
    return runner
