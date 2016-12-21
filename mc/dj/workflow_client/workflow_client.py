import requests


class MissionControlWorkflowClient(object):
    def __init__(self, base_url=None, request_client=None):
        self.base_url = base_url
        self.request_client = request_client or requests

    def fetch_workflows(self, query_params=None):
        workflows_url = self.base_url + 'workflows/'
        if query_params:
            response = self.request_client.get(workflows_url, query_params)
        else:
            response = self.request_client.get(workflows_url)
        return response.json()

    def claim_workflows(self, uuids=None):
        claim_workflows_url = self.base_url + 'claim_workflows/'
        response = self.request_client.post(claim_workflows_url, {'uuids': uuids})
        return response.json()

    def update_workflows(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            workflow_url = self.base_url + 'workflows/' + _uuid + '/'
            response = self.request_client.patch(workflow_url, updates_for_uuid)
            results_by_uuid[_uuid] = response.json()
        return results_by_uuid
