import requests
import logging


class MissionControlFlowClient(object):
    def __init__(self, base_url=None, request_client=None, logger=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests
        self.logger = logger or logging

    def generate_urls(self):
        return {
            'flows': self.base_url + 'flows/',
            'claim_flows': self.base_url + 'claim_flows/'
        }

    def create_flow(self, flow=None):
        response = self.request_client.post(self.urls['flows'], json=flow)
        return self.json_raise_for_status(response=response)

    def json_raise_for_status(self, response=None):
        try:
            response.raise_for_status()
            return response.json()
        except Exception as exception:
            wrapped_exception = Exception("%s; %s" % (exception,
                                                      response.content))
            self.logger.exception(wrapped_exception)
            raise wrapped_exception

    def fetch_flows(self, query_params=None):
        args = [self.urls['flows']]
        if query_params: args.append(query_params)
        response = self.request_client.get(*args)
        return self.json_raise_for_status(response=response)

    def fetch_flow_by_uuid(self, uuid=None):
        fetch_flows_result = self.fetch_flows(query_params={'uuid': uuid})
        if len(fetch_flows_result) == 1: return fetch_flows_result[0]
        else: return None

    def fetch_tickable_flows(self):
        return self.fetch_flows(query_params={'claimed': False,
                                              'tickable': True})

    def claim_flows(self, uuids=None):
        response = self.request_client.post(self.urls['claim_flows'], 
                                            json={'uuids': uuids})
        return self.json_raise_for_status(response=response)

    def update_flows(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            flow_url = self.base_url + 'flows/' + _uuid + '/'
            response = self.request_client.patch(flow_url,
                                                 json=updates_for_uuid)
            try:
                result = self.json_raise_for_status(response=response)
            except Exception as exception:
                result = {"error", str(exception)}
            results_by_uuid[_uuid] = result
        return results_by_uuid
