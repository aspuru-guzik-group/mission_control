import requests


class MissionControlFlowClient(object):
    def __init__(self, base_url=None, request_client=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests

    def generate_urls(self):
        return {
            'flows': self.base_url + 'flows/',
            'claim_flows': self.base_url + 'claim_flows/'
        }

    def create_flow(self, flow=None):
        response = self.request_client.post(self.urls['flows'], data=flow)
        return response.json()

    def fetch_flows(self, query_params=None):
        if query_params:
            response = self.request_client.get(self.urls['flows'], query_params)
        else:
            response = self.request_client.get(self.urls['flows'])
        return response.json()

    def fetch_flow_by_uuid(self, uuid=None):
        fetch_flows_result = self.fetch_flows(query_params={'uuid': uuid})
        return fetch_flows_result.get(uuid, None)

    def claim_flows(self, uuids=None):
        response = self.request_client.post(self.urls['claim_flows'], 
                                            {'uuids': uuids})
        return response.json()

    def update_flows(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            flow_url = self.base_url + 'flows/' + _uuid + '/'
            response = self.request_client.patch(flow_url, updates_for_uuid)
            results_by_uuid[_uuid] = response.json()
        return results_by_uuid
