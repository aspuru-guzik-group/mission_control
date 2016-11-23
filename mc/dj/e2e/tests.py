import copy
from django.conf.urls import url, include
from django.test import TestCase, override_settings
import json
from job_runners.base_job_runner import BaseJobRunner
from jobs.models import Job
from job_spec_client.job_spec_client import MissionControlJobSpecClient

BASE_PATH = 'test_api'

urlpatterns = [
    url(r'^%s' % BASE_PATH, include('jobs.urls')),
]

@override_settings(ROOT_URLCONF=__name__)
class ClientJobRunner_X_Api_TestCase(TestCase):
    def setUp(self):
        self.num_jobs = 10
        self.jobs = self.generate_jobs()
        self.runner = self.generate_runner()
        orig_patch = self.client.patch
        def json_patch(path, data=None, **kwargs):
            return orig_patch(path, json.dumps(data),
                              content_type='application/json', **kwargs)
        self.client.patch = json_patch

    def generate_jobs(self):
        jobs = []
        for i in range(self.num_jobs):
            jobs.append(Job.objects.create())

    def generate_runner(self):
        class StubJobDirFactory(object):
            def build_dir_for_spec(self, job_spec=None):
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
            job_spec_client=MissionControlJobSpecClient(
                base_url='/%s' % BASE_PATH,
                request_client=self.client
            ),
            job_dir_factory=StubJobDirFactory(),
            transfer_client=StubTransferClient()
        )
        return runner

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
            if j.status == self.runner.job_spec_client.Statuses.Claimed.name}
        state['pending_jobs'] = {
            key:j for key, j in state['jobs'].items()
            if j.status == self.runner.job_spec_client.Statuses.Pending.name}
        state['completed_jobs'] = {
            key:j for key, j in state['jobs'].items()
            if j.status == self.runner.job_spec_client.Statuses.Completed.name}
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
