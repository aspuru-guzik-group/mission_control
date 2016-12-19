import json
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..models import Workflow, WorkflowStatuses
from ..serializers import WorkflowSerializer


@override_settings(ROOT_URLCONF='workflows.urls')
class ListWorkflowsTestCase(APITestCase):
    def test_list_workflows(self):
        workflows = [Workflow.objects.create(name="workflow_%s" % i)
                     for i in range(3)]
        response = self.client.get('/workflows/')
        expected_data = [WorkflowSerializer(workflow).data
                         for workflow in workflows]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

    def test_status_filtering(self):
        statuses = [WorkflowStatuses.Pending.name,
                    WorkflowStatuses.Claimed.name]
        workflows_by_status = {
            status: [Workflow.objects.create(name="workflow_%s" % i,
                                             status=status)
                     for i in range(3)]
            for status in statuses
        }
        for status in statuses:
            response = self.client.get('/workflows/', {'status': status})
            expected_data = [WorkflowSerializer(workflow).data
                             for workflow in workflows_by_status[status]]
            self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                             sorted(expected_data, key=lambda j: j['uuid']))

@override_settings(ROOT_URLCONF='workflows.urls')
class PatchWorkflowTestCase(APITestCase):
    def setUp(self):
        self.workflows = [Workflow.objects.create(name="workflow_%s" % i)
                     for i in range(1)]

    def test_patch_workflow(self):
        workflow_to_patch = self.workflows[0]
        new_values = {'name': 'new_name'}
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

@override_settings(ROOT_URLCONF='workflows.urls')
class ClaimWorkflowTestCase(TestCase):
    def setUp(self):
        self.unclaimed_workflows = [
            Workflow.objects.create(name="workflow_%s" % i) for i in range(1)]
        self.claimed_workflows = [Workflow.objects.create(
            name="workflow_%s" % i, status=WorkflowStatuses.Claimed.name)
            for i in range(1)]
        self.all_workflows = self.unclaimed_workflows + self.claimed_workflows
        self.workflows_to_claim = [j for j in self.unclaimed_workflows[:-1]] + [
            j for j in self.claimed_workflows[:-1]]
        self._claimWorkflows()

    def _claimWorkflows(self):
        csv_uuids = ','.join([workflow.uuid
                              for workflow in self.workflows_to_claim])
        self.response = self.client.post('/claim_workflows/',
                                         {'uuids': csv_uuids})

    def test_response_data(self):
        self.assertEqual(self.response.status_code, 200)
        expected_data = {str(workflow.uuid): WorkflowSerializer(workflow).data
                         for workflow in self.workflows_to_claim}
        self.assertEqual(json.loads(self.response.content.decode()),
                         expected_data)

    def test_db_state(self):
        statuses = {workflow.uuid: workflow.status
                    for workflow in Workflow.objects.all()}
        expected_statuses = {}
        for workflow in self.all_workflows:
            combined = set(self.workflows_to_claim + self.claimed_workflows)
            if workflow in combined:
                expected_status = WorkflowStatuses.Claimed.name
            else:
                expected_status = WorkflowStatuses.Pending.name
            expected_statuses[str(workflow.uuid)] = expected_status
        self.assertEqual(statuses, expected_statuses)

