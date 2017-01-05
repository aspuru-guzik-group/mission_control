import requests


class MissionControlFlowClient(object):
    def __init__(self, base_url=None, request_client=None):
        self.base_url = base_url
        self.request_client = request_client or requests

    def fetch_flows(self, query_params=None):
        flows_url = self.base_url + 'flows/'
        if query_params:
            response = self.request_client.get(flows_url, query_params)
        else:
            response = self.request_client.get(flows_url)
        return response.json()

    def claim_flows(self, uuids=None):
        claim_flows_url = self.base_url + 'claim_flows/'
        response = self.request_client.post(claim_flows_url, {'uuids': uuids})
        return response.json()

    def update_flows(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            flow_url = self.base_url + 'flows/' + _uuid + '/'
            response = self.request_client.patch(flow_url, updates_for_uuid)
            results_by_uuid[_uuid] = response.json()
        return results_by_uuid
