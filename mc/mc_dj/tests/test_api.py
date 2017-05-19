import json
from unittest.mock import call, patch

from django.conf.urls import url, include
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..constants import JobStatuses
from .. import models as _models
from .. import serializers as _serializers
from .. import urls as _urls
from .. import views as _views


def key_by_uuid(items=None):
    return {item['uuid']: item for item in items}

BASE_PATH = 'test_api/'
BASE_URL = '/' + BASE_PATH
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_urls.__name__))
]

@override_settings(ROOT_URLCONF=__name__)
class BaseTestCase(APITestCase): pass

class ListFlowsTestCase(BaseTestCase):
    def test_list_flows(self):
        flows = [_models.Flow.objects.create() for i in range(3)]
        response = self.client.get(BASE_URL + 'flows/')
        expected_data = [_serializers.FlowSerializer(flow).data
                         for flow in flows]
        self.assertEqual(key_by_uuid(response.data), key_by_uuid(expected_data))

    def test_status_filtering(self):
        statuses = [_models.FlowStatuses.PENDING.name,
                    _models.FlowStatuses.RUNNING.name]
        flows_by_status = {
            status: [_models.Flow.objects.create(status=status)
                     for i in range(3)]
            for status in statuses
        }
        for status in statuses:
            response = self.client.get(BASE_URL + 'flows/', {'status': status})
            expected_data = [_serializers.FlowSerializer(flow).data
                             for flow in flows_by_status[status]]
            self.assertEqual(key_by_uuid(response.data),
                             key_by_uuid(expected_data))

    def test_uuid_filtering(self):
        flows = [_models.Flow.objects.create() for i in range(2)]
        flows_by_uuid = {flow.uuid: flow for flow in flows}
        for uuid, flow in flows_by_uuid.items():
            response = self.client.get(BASE_URL + 'flows/', {'uuid': uuid})
            expected_data = [_serializers.FlowSerializer(flow).data]
            self.assertEqual(key_by_uuid(response.data),
                             key_by_uuid(expected_data))

class TickableFlowsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_groups = self.generate_flow_groups()

    def generate_flow_groups(self):
        flow_groups = {
            'tickable': [
                _models.Flow.objects.create(claimed=True, status=status.name)
                for status in _models.FlowStatuses.tickable_statuses
            ],
            'untickable': [
                _models.Flow.objects.create(
                    claimed=True, status=_models.FlowStatuses.COMPLETED.name)
            ],
        }
        for flow_group in flow_groups.items():
            self.assertTrue(len(flow_group) > 0)
        return flow_groups

    def test_returns_tickable_flows(self):
        response = self.client.get(BASE_URL + 'flows/', {'tickable': True})
        expected_data = [
            _serializers.FlowSerializer(flow).data
            for flow in self.flow_groups['tickable']
        ]
        self.assertEqual(key_by_uuid(response.data), key_by_uuid(expected_data))

class PatchFlowTestCase(BaseTestCase):
    def setUp(self):
        self.flows = [_models.Flow.objects.create() for i in range(1)]

    def test_patch_flow(self):
        flow_to_patch = self.flows[0]
        new_values = {'status': _models.FlowStatuses.RUNNING.name}
        response = self.client.patch(
            BASE_URL + 'flows/%s/' % flow_to_patch.uuid, new_values)
        self.assertEqual(response.status_code, 200)
        patched_flow = _models.Flow.objects.get(uuid=flow_to_patch.uuid)
        patched_flow_attrs = {attr: getattr(patched_flow, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_flow_attrs, new_values)

    def test_patch_with_bad_status(self):
        response = self.client.patch(
            BASE_URL + 'flows/%s/' % self.flows[0].uuid, {'status': 'BADDO!'})
        self.assertEqual(response.status_code, 400)

class ClaimFlowTestCase(BaseTestCase):
    def setUp(self):
        self.unclaimed_flows = [_models.Flow.objects.create() for i in range(1)]
        self.claimed_flows = [_models.Flow.objects.create(claimed=True)
                                  for i in range(1)]
        self.all_flows = self.unclaimed_flows + self.claimed_flows
        self.flows_to_claim = [j for j in self.unclaimed_flows[:-1]] + \
                [j for j in self.claimed_flows[:-1]]
        self._claimFlows()

    def _claimFlows(self):
        uuids = [flow.uuid for flow in self.flows_to_claim]
        self.response = self.client.post(
            BASE_URL + 'claim_flows/', 
            json.dumps({'uuids': uuids}),
            content_type='application/json'
        )

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {flow.uuid: _serializers.FlowSerializer(flow).data
                         for flow in self.flows_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        claimed_values = {flow.uuid: flow.claimed
                          for flow in _models.Flow.objects.all()}
        expected_claimed_values = {}
        for flow in self.all_flows:
            combined = set(self.flows_to_claim + self.claimed_flows)
            if flow in combined:
                expected_claimed_value = True
            else:
                expected_claimed_value = False
            expected_claimed_values[flow.uuid] = expected_claimed_value
        self.assertEqual(claimed_values, expected_claimed_values)


class ListJobsTestCase(BaseTestCase):
    def test_list_jobs(self):
        jobs = [_models.Job.objects.create(label="job_%s" % i)
                for i in range(3)]
        response = self.client.get(BASE_URL + 'jobs/')
        expected_data = [_serializers.JobSerializer(job).data for job in jobs]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

    def test_status_filtering(self):
        statuses = [JobStatuses.PENDING.name, JobStatuses.RUNNING.name]
        jobs_by_status = {
            status: [
                _models.Job.objects.create(label="job_%s" % i, status=status)
                for i in range(3)
            ]
            for status in statuses
        }
        for status in statuses:
            response = self.client.get(BASE_URL + 'jobs/', {'status': status})
            expected_data = [_serializers.JobSerializer(job).data
                             for job in jobs_by_status[status]]
            self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                             sorted(expected_data, key=lambda j: j['uuid']))

    def test_uuid_filtering(self):
        jobs = [_models.Job.objects.create() for i in range(2)]
        jobs_by_uuid = {job.uuid: job for job in jobs}
        for uuid, job in jobs_by_uuid.items():
            response = self.client.get(BASE_URL + 'jobs/', {'uuid': uuid})
            expected_data = [_serializers.JobSerializer(job).data]
            self.assertEqual(self.key_by_uuid(response.data),
                             self.key_by_uuid(expected_data))

    def key_by_uuid(self, items=None):
        return {item['uuid']: item for item in items}

class PatchJobTestCase(BaseTestCase):
    def setUp(self):
        self.jobs = [_models.Job.objects.create(label="job_%s" % i)
                     for i in range(1)]

    def test_patch_job(self):
        job_to_patch = self.jobs[0]
        new_values = {'label': 'new_label', 'data': json.dumps({'new': 'data'})}
        response = self.client.patch(BASE_URL + 'jobs/%s/' % job_to_patch.uuid,
                                     new_values)
        self.assertEqual(response.status_code, 200)
        patched_job = _models.Job.objects.get(uuid=job_to_patch.uuid)
        patched_job_attrs = {attr: getattr(patched_job, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_job_attrs, new_values)

    def test_patch_with_bad_status(self):
        response = self.client.patch(BASE_URL + 'jobs/%s/' % self.jobs[0].uuid,
                                     {'status': 'BADDO!'})
        self.assertEqual(response.status_code, 400)

class ClaimJobTestCase(BaseTestCase):
    def setUp(self):
        self.unclaimed_jobs = [_models.Job.objects.create(label="job_%s" % i)
                               for i in range(1)]
        self.claimed_jobs = [_models.Job.objects.create(
            label="job_%s" % i, status=JobStatuses.RUNNING.name)
            for i in range(1)]
        self.all_jobs = self.unclaimed_jobs + self.claimed_jobs
        self.jobs_to_claim = [j for j in self.unclaimed_jobs[:-1]] + [
            j for j in self.claimed_jobs[:-1]]
        self._claimJobs()

    def _claimJobs(self):
        uuids = [job.uuid for job in self.jobs_to_claim]
        self.response = self.client.post(BASE_URL + 'claim_jobs/',
                                         json.dumps({'uuids': uuids}),
                                         content_type='application/json')

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {str(job.uuid): _serializers.JobSerializer(job).data
                         for job in self.jobs_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        statuses = {job.uuid: job.status for job in _models.Job.objects.all()}
        expected_statuses = {}
        for job in self.all_jobs:
            if job in set(self.jobs_to_claim + self.claimed_jobs):
                expected_status = JobStatuses.RUNNING.name
            else:
                expected_status = JobStatuses.PENDING.name
            expected_statuses[str(job.uuid)] = expected_status
        self.assertEqual(statuses, expected_statuses)

class FlushTestCase(BaseTestCase):
    def setUp(self):
        self.create_models()

    def create_models(self):
        for model in _models.missions_models:
            for i in range(3): model.objects.create()

    def _flush(self):
        return self.client.get(BASE_URL + 'flush/')

    def test_response_data(self):
        response = self._flush()
        self.assertEqual(response.status_code, 200)
        for model in _models.missions_models:
            self.assertEqual(model.objects.count(), 0)
        expected_flush_results = {model.__name__: 'flushed'
                                  for model in _models.missions_models}
        self.assertEqual(json.loads(response.content.decode()), 
                         expected_flush_results)

class QueueItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue = _models.Queue.objects.create()
        self.mock_serialized_items = ['item_%s' % i for i in range(3)]
        self.expected_models = {}
        self.expected_serializers = {}
        for item_type in ['Job', 'Flow']:
            self.expected_models[item_type] = getattr(_models, item_type)
            self.expected_serializers[item_type] = getattr(
                _serializers, '%sSerializer' % item_type)

    @patch.object(_views, '_queue_utils')
    def test_returns_serialized_claimed_queue_items(self, mock_queue_utils):
        mock_queue_utils.serialize_queue_items.return_value = \
                self.mock_serialized_items
        url = '{base_url}queues/{queue_key}/claim_items/'.format(
            base_url=BASE_URL, queue_key=self.queue.uuid)
        response = self.client.post(url)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(mock_queue_utils.claim_queue_items.call_args,
                         call(queue=self.queue, models=self.expected_models))
        self.assertEqual(
            mock_queue_utils.serialize_queue_items.call_args,
            call(queue=self.queue,
                 items=mock_queue_utils.claim_queue_items.return_value,
                 serializers=self.expected_serializers)
        )
        self.assertEqual(result['items'], self.mock_serialized_items)
