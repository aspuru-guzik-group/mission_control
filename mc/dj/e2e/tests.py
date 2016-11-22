from django.conf.urls import url, include
from django.test import TestCase, override_settings
from local_job_runner.base_daemon import BaseDaemon
from jobs.models import Job
from job_spec_client import MissionControlJobSpecClient

BASE_URL = 'test_api'

urlpatterns = [
    url(r'^%s' % BASE_URL, include('jobs.urls')),
]

@override_settings(ROOT_URLCONF='__main__')
class ClientDaemon_X_Api_TestCase(TestCase):
    def setUp(self):
        self.num_jobs = 10
        self.jobs = self.generate_jobs()
        self.runner = self.generate_runner()

    def generate_jobs(self):
        jobs = []
        for i in range(self.num_jobs):
            jobs.append(Job.objects.create())

    def generate_runner(self):
        class StubTransferClient(object):
            pass

        class StubJobDirFactory(object):
            def build_dir_for_spec(self, job_spec=None):
                dir_meta = {'dir': 'test_dir', 'entrypoint': 'job.sh'}
                return dir_meta

        runner = BaseDaemon(
            job_spec_client=MissionControlJobSpecClient(
                base_url=BASE_URL + '/'
            ),
            job_dir_factory=StubJobDirFactory(),
            transfer_client=StubTransferClient()
        )
        return runner

    def test_job_cycle(self):
        # Tick 1: should request and claim jobs.
        # Tick 2: should be executing jobs, should not be fetching new jobs
        # (running jobs maxed out)
        # <set jobs to finish executing>
        # Tick 3: should start transferring jobs, fetching new jobs.
        # <set jobs to finish transferring>
        # Tick 4: should update job status.
        pass
        self.fail()


