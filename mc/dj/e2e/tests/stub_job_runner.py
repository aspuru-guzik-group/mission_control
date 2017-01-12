from job_runners.base_job_runner import BaseJobRunner
from job_client.job_client import MissionControlJobClient


def generate_stub_job_runner(base_url=None, request_client=None):
    class StubJobDirFactory(object):
        def build_dir_for_spec(self, job=None):
            return {}

    class StubTransferClient(object):
        def start_transfer(self, job=None):
            return {}

        def get_transfer_state(self, job=None):
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
        job_dir_factory=StubJobDirFactory(),
        transfer_client=StubTransferClient()
    )
    return runner
