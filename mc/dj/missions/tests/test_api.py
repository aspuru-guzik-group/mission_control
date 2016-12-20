import json
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..models import Workflow, WorkflowStatuses
from ..serializers import WorkflowSerializer


@override_settings(ROOT_URLCONF='missions.urls')
class ListWorkflowsTestCase(APITestCase):
    def test_list_workflows(self):
        workflows = [Workflow.objects.create() for i in range(3)]
        response = self.client.get('/workflows/')
        expected_data = [WorkflowSerializer(workflow).data
                         for workflow in workflows]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

    def test_status_filtering(self):
        statuses = [WorkflowStatuses.PENDING.name,
                    WorkflowStatuses.RUNNING.name]
        workflows_by_status = {
            status: [Workflow.objects.create(status=status) for i in range(3)]
            for status in statuses
        }
        for status in statuses:
            response = self.client.get('/workflows/', {'status': status})
            expected_data = [WorkflowSerializer(workflow).data
                             for workflow in workflows_by_status[status]]
            self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                             sorted(expected_data, key=lambda j: j['uuid']))

@override_settings(ROOT_URLCONF='missions.urls')
class PatchWorkflowTestCase(APITestCase):
    def setUp(self):
        self.workflows = [Workflow.objects.create() for i in range(1)]

    def test_patch_workflow(self):
        workflow_to_patch = self.workflows[0]
        new_values = {'status': WorkflowStatuses.RUNNING.name}
        response = self.client.patch('/workflows/%s/' % workflow_to_patch.uuid,
                                     new_values)
        self.assertEqual(response.status_code, 200)
        patched_workflow = Workflow.objects.get(uuid=workflow_to_patch.uuid)
        patched_workflow_attrs = {attr: getattr(patched_workflow, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_workflow_attrs, new_values)

    def test_patch_with_bad_status(self):
        response = self.client.patch('/workflows/%s/' % self.workflows[0].uuid,
                                     {'status': 'BADDO!'})
        self.assertEqual(response.status_code, 400)

@override_settings(ROOT_URLCONF='missions.urls')
class ClaimWorkflowTestCase(TestCase):
    def setUp(self):
        self.unclaimed_workflows = [Workflow.objects.create() for i in range(1)]
        self.claimed_workflows = [Workflow.objects.create(claimed=True)
                                  for i in range(1)]
        self.all_workflows = self.unclaimed_workflows + self.claimed_workflows
        self.workflows_to_claim = [j for j in self.unclaimed_workflows[:-1]] + \
                [j for j in self.claimed_workflows[:-1]]
        self._claimWorkflows()

    def _claimWorkflows(self):
        csv_uuids = ','.join([workflow.uuid
                              for workflow in self.workflows_to_claim])
        self.response = self.client.post('/claim_workflows/',
                                         {'uuids': csv_uuids})

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {workflow.uuid: WorkflowSerializer(workflow).data
                         for workflow in self.workflows_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        claimed_values = {workflow.uuid: workflow.claimed
                          for workflow in Workflow.objects.all()}
        expected_claimed_values = {}
        for workflow in self.all_workflows:
            combined = set(self.workflows_to_claim + self.claimed_workflows)
            if workflow in combined:
                expected_claimed_value = True
            else:
                expected_claimed_value = False
            expected_claimed_values[workflow.uuid] = expected_claimed_value
        self.assertEqual(claimed_values, expected_claimed_values)

