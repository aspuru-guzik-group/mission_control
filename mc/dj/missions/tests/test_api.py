import json
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..models import Flow, FlowStatuses
from ..serializers import FlowSerializer


def key_by_uuid(items=None):
    return {item['uuid']: item for item in items}

@override_settings(ROOT_URLCONF='missions.urls')
class ListFlowsTestCase(APITestCase):
    def test_list_flows(self):
        flows = [Flow.objects.create() for i in range(3)]
        response = self.client.get('/flows/')
        expected_data = [FlowSerializer(flow).data
                         for flow in flows]
        self.assertEqual(key_by_uuid(response.data), key_by_uuid(expected_data))

    def test_status_filtering(self):
        statuses = [FlowStatuses.PENDING.name,
                    FlowStatuses.RUNNING.name]
        flows_by_status = {
            status: [Flow.objects.create(status=status) for i in range(3)]
            for status in statuses
        }
        for status in statuses:
            response = self.client.get('/flows/', {'status': status})
            expected_data = [FlowSerializer(flow).data
                             for flow in flows_by_status[status]]
            self.assertEqual(key_by_uuid(response.data),
                             key_by_uuid(expected_data))

@override_settings(ROOT_URLCONF='missions.urls')
class TickableFlowsTestCase(APITestCase):
    def setUp(self):
        super().setUp()
        self.flow_groups = self.generate_flow_groups()

    def generate_flow_groups(self):
        flow_groups = {
            'tickable': [
                Flow.objects.create(claimed=True, status=status.name)
                for status in FlowStatuses.tickable_statuses
            ],
            'untickable': [
                Flow.objects.create(claimed=True,
                                        status=FlowStatuses.COMPLETED.name)
            ],
        }
        for flow_group in flow_groups.items():
            self.assertTrue(len(flow_group) > 0)
        return flow_groups

    def test_returns_tickable_flows(self):
        response = self.client.get('/flows/', {'tickable': True})
        expected_data = [
            FlowSerializer(flow).data
            for flow in self.flow_groups['tickable']
        ]
        self.assertEqual(key_by_uuid(response.data), key_by_uuid(expected_data))

@override_settings(ROOT_URLCONF='missions.urls')
class PatchFlowTestCase(APITestCase):
    def setUp(self):
        self.flows = [Flow.objects.create() for i in range(1)]

    def test_patch_flow(self):
        flow_to_patch = self.flows[0]
        new_values = {'status': FlowStatuses.RUNNING.name}
        response = self.client.patch('/flows/%s/' % flow_to_patch.uuid,
                                     new_values)
        self.assertEqual(response.status_code, 200)
        patched_flow = Flow.objects.get(uuid=flow_to_patch.uuid)
        patched_flow_attrs = {attr: getattr(patched_flow, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_flow_attrs, new_values)

    def test_patch_with_bad_status(self):
        response = self.client.patch('/flows/%s/' % self.flows[0].uuid,
                                     {'status': 'BADDO!'})
        self.assertEqual(response.status_code, 400)

@override_settings(ROOT_URLCONF='missions.urls')
class ClaimFlowTestCase(TestCase):
    def setUp(self):
        self.unclaimed_flows = [Flow.objects.create() for i in range(1)]
        self.claimed_flows = [Flow.objects.create(claimed=True)
                                  for i in range(1)]
        self.all_flows = self.unclaimed_flows + self.claimed_flows
        self.flows_to_claim = [j for j in self.unclaimed_flows[:-1]] + \
                [j for j in self.claimed_flows[:-1]]
        self._claimFlows()

    def _claimFlows(self):
        csv_uuids = ','.join([flow.uuid
                              for flow in self.flows_to_claim])
        self.response = self.client.post('/claim_flows/',
                                         {'uuids': csv_uuids})

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {flow.uuid: FlowSerializer(flow).data
                         for flow in self.flows_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        claimed_values = {flow.uuid: flow.claimed
                          for flow in Flow.objects.all()}
        expected_claimed_values = {}
        for flow in self.all_flows:
            combined = set(self.flows_to_claim + self.claimed_flows)
            if flow in combined:
                expected_claimed_value = True
            else:
                expected_claimed_value = False
            expected_claimed_values[flow.uuid] = expected_claimed_value
        self.assertEqual(claimed_values, expected_claimed_values)

