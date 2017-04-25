import json
import requests
import logging


class MissionControlClient(object):
    def __init__(self, base_url=None, request_client=None, logger=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests
        self.logger = logger or logging

    def generate_urls(self):
        return {
            'flows': self.base_url + 'flows/',
            'claim_flows': self.base_url + 'claim_flows/',
            'jobs': self.base_url + 'jobs/',
            'claim_jobs': self.base_url + 'claim_jobs/',
            'flush': self.base_url + 'flush/',
            'queues': self.base_url + 'queues/',
        }

    def create_flow(self, flow_kwargs=None):
        formatted_flow_kwargs = self.format_flow_kwargs(flow_kwargs)
        response = self.request_client.post(self.urls['flows'],
                                            json=formatted_flow_kwargs)
        return self.format_fetched_flow(
            fetched_flow=self.json_raise_for_status(response=response))

    def format_flow_kwargs(self, flow_kwargs=None):
        formatted_flow_kwargs = {**flow_kwargs}
        if 'serialization' in flow_kwargs:
            formatted_flow_kwargs['serialization'] = json.dumps(
                flow_kwargs['serialization'])
        return formatted_flow_kwargs

    def format_fetched_flow(self, fetched_flow=None):
        formatted_flow = {
            **fetched_flow,
            'serialization': json.loads(fetched_flow.get('serialization', '{}'))
        }
        return formatted_flow

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
        return self.format_flow_list_response(response=response)

    def format_flow_list_response(self, response=None):
        return [
            self.format_fetched_flow(fetched_flow)
            for fetched_flow in self.json_raise_for_status(response=response)
        ]

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
        claimed_flows_by_uuid = self.json_raise_for_status(response=response)
        return {
            uuid: self.format_fetched_flow(claimed_flow)
            for uuid, claimed_flow in claimed_flows_by_uuid.items()
        }

    def update_flows(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            formatted_updates_for_uuid = self.format_flow_kwargs(
                flow_kwargs=updates_for_uuid)
            flow_url = self.base_url + 'flows/' + _uuid + '/'
            response = self.request_client.patch(
                flow_url, json=formatted_updates_for_uuid)
            try:
                result = self.json_raise_for_status(response=response)
            except Exception as exception:
                result = {"error", str(exception)}
            results_by_uuid[_uuid] = result
        return results_by_uuid

    def create_job(self, job_kwargs=None):
        formatted_job_kwargs = self.format_job_kwargs(job_kwargs)
        response = self.request_client.post(self.urls['jobs'],
                                            json=formatted_job_kwargs)
        return self.format_fetched_job(
            fetched_job=self.json_raise_for_status(response=response))

    def format_job_kwargs(self, job_kwargs=None):
        formatted_job_kwargs = {**job_kwargs}
        if 'data' in job_kwargs:
            formatted_job_kwargs['data'] = self.serialize_job_data(
                job_kwargs['data'])
        return formatted_job_kwargs

    def serialize_job_data(self, job_data=None):
        return json.dumps(job_data)

    def fetch_jobs(self, query_params=None):
        args = [self.urls['jobs']]
        if query_params: args.append(query_params)
        response = self.request_client.get(*args)
        jobs = [
            self.format_fetched_job(fetched_job)
            for fetched_job in self.json_raise_for_status(response=response)
        ]
        return jobs

    def format_fetched_job(self, fetched_job=None):
        formatted_job = {
            **fetched_job,
            'data': self.deserialize_job_data(fetched_job.get('data', '{}'))
        }
        return formatted_job

    def deserialize_job_data(self, serialized_job_data=None):
        return json.loads(serialized_job_data)

    def fetch_claimable_jobs(self):
        return self.fetch_jobs(query_params={'status': 'PENDING'})

    def fetch_job_by_uuid(self, uuid=None):
        result = self.fetch_jobs(query_params={'uuid': uuid})
        if len(result) == 1: return result[0]
        else: raise Exception("Could not get job with uuid '%s'" % uuid)

    def claim_jobs(self, uuids=None):
        response = self.request_client.post(self.urls['claim_jobs'],
                                            json={'uuids': uuids})
        claim_results = self.json_raise_for_status(response=response)
        formatted_claim_results = {}
        for uuid, claimed_job in claim_results.items():
            if claimed_job: claimed_job = self.format_fetched_job(claimed_job)
            formatted_claim_results[uuid] = claimed_job
        return formatted_claim_results

    def update_jobs(self, updates_by_uuid=None):
        results_by_uuid = {}
        for uuid, updates_for_uuid in updates_by_uuid.items():
            results_by_uuid[uuid] = self.update_job(uuid=uuid,
                                                    updates=updates_for_uuid)
        return results_by_uuid

    def update_job(self, uuid=None, updates=None):
        job_url = self.urls['jobs'] + uuid + '/'
        formatted_job_kwargs = self.format_job_kwargs(job_kwargs=updates)
        response = self.request_client.patch(job_url,
                                             json=formatted_job_kwargs)
        return self.json_raise_for_status(response=response)

    def flush_mc_db(self):
        url = self.urls['flush']
        response = self.request_client.get(url)
        return self.json_raise_for_status(response=response)

    def create_queue(self, queue_kwargs=None):
        response = self.request_client.post(self.urls['queues'],
                                            json=queue_kwargs)
        return self.json_raise_for_status(response=response)

    def claim_queue_items(self, queue_key=None, params=None,
                          item_transform_fn=None):
        params = params or {}
        url = '{queues_root}{queue_key}/claim_items/'.format(
            queues_root=self.urls['queues'], queue_key=queue_key)
        response = self.request_client.post(url, data=json.dumps(params))
        result = self.json_raise_for_status(response=response)
        if item_transform_fn:
            result['items'] = [item_transform_fn(item)
                               for item in result['items']]
        return result

    def claim_job_queue_items(self, *args, **kwargs):
        return self.claim_queue_items(
            *args, item_transform_fn=self.format_fetched_job, **kwargs)

    def claim_flow_queue_items(self, *args, **kwargs):
        return self.claim_queue_items(
            *args, item_transform_fn=self.format_fetched_flow, **kwargs)
