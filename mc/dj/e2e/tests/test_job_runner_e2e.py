import copy
from django.conf.urls import url, include
from django.test import TestCase, override_settings
from jobs.models import Job

from mc_utils import test_utils
from .stub_job_runner import generate_stub_job_runner

BASE_PATH = 'test_api/'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('jobs.urls')),
]

@override_settings(ROOT_URLCONF=__name__)
class JobRunnerE2ETestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client(request_client=self.client)
        self.num_jobs = 10
        self.jobs = self.generate_jobs()
        self.runner = generate_stub_job_runner(base_url='/%s' % BASE_PATH,
                                               request_client=self.client)

    def generate_jobs(self):
        jobs = []
        for i in range(self.num_jobs):
            jobs.append(Job.objects.create())

    def test_job_cycle(self):
        states = {0: self.get_current_state()}
        for i in range(1, 10):
            self.runner.tick()
            prev_state = states[i - 1]
            states[i] = self.get_current_state(prev_state=prev_state)
            self.check_state(state=states[i], prev_state=prev_state)

    def get_current_state(self, prev_state=None):
        state = {'tick_counter': self.runner.tick_counter}
        state['jobs'] = {j.uuid: j for j in Job.objects.all()}
        state['claimed_jobs'] = {
            key:j for key, j in state['jobs'].items()
            if j.status == self.runner.job_client.Statuses.RUNNING.name}
        state['pending_jobs'] = {
            key:j for key, j in state['jobs'].items()
            if j.status == self.runner.job_client.Statuses.PENDING.name}
        state['completed_jobs'] = {
            key:j for key, j in state['jobs'].items()
            if j.status == self.runner.job_client.Statuses.COMPLETED.name}
        state['executing_jobs'] = copy.deepcopy(self.runner.executing_jobs)
        state['transferring_jobs'] = copy.deepcopy(
            self.runner.transferring_jobs)
        if prev_state:
            state['newly_claimed_jobs'] = {
                key:j for key, j in state['claimed_jobs'].items()
                if key not in prev_state['claimed_jobs']
            }
            state['newly_completed_jobs'] = {
                key:j for key, j in state['completed_jobs'].items()
                if key not in prev_state['completed_jobs']
            }
            state['newly_executing_jobs'] = {
                key:j for key, j in state['executing_jobs'].items()
                if key not in prev_state['executing_jobs']
            }
            state['newly_executed_jobs'] = {
                key:j for key, j in prev_state['executing_jobs'].items()
                if key not in state['executing_jobs']
            }
            state['newly_transferring_jobs'] = {
                key:j for key, j in state['transferring_jobs'].items()
                if key not in prev_state['transferring_jobs']
            }
            state['newly_transferred_jobs'] = {
                key:j for key, j in prev_state['transferring_jobs'].items()
                if key not in state['transferring_jobs']
            }
        return state

    def check_state(self, state=None, prev_state=None):
        self.assertEqual(
            len(state.get('newly_transferring_jobs', {})),
            len(state.get('newly_executed_jobs', {}))
        )
        self.assertEqual(
            len(state.get('newly_executing_jobs', {})),
            len(state.get('newly_claimed_jobs', {}))
        )
        if state['tick_counter'] == 1:
            expected_execution_slots = self.runner.max_executing_jobs
        else:
            expected_execution_slots = \
                    len(state.get('newly_executed_jobs', {}))
        self.assertEqual(len(state.get('newly_claimed_jobs', {})),
                         min(len(prev_state['pending_jobs']),
                             expected_execution_slots))
        self.assertEqual(
            len(state.get('newly_completed_jobs', {})),
            len(state.get('newly_transferred_jobs', {}))
        )
