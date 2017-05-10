import json
import requests as _requests
import logging


class HttpDao(object):
    def __init__(self, base_url=None, requests=None, logger=None,
                 max_error_len=1024):
        self.base_url = base_url
        self.requests = requests or _requests
        self.logger = logger or logging
        self.max_error_len = max_error_len

        self.urls = self.generate_urls()

    def generate_urls(self):
        return {
            'flows': self.base_url + 'flows/',
            'flush': self.base_url + 'flush/',
            'jobs': self.base_url + 'jobs/',
            'queues': self.base_url + 'queues/',
        }

    def create_flow(self, flow=None):
        response = self.requests.post(self.urls['flows'], json=flow)
        return self.json_raise_for_status(response=response)

    def json_raise_for_status(self, response=None):
        try:
            response.raise_for_status()
            return response.json()
        except Exception as exception:
            wrapped_exception = Exception("%s; %s" % (
                exception, self._elide_string(response.content.decode(),
                                              max_len=self.max_error_len)))
            self.logger.exception(wrapped_exception)
            raise wrapped_exception

    def _elide_string(self, string=None, max_len=None):
        if max_len is None: return string
        if len(string) > max_len: return string[:max_len] + '...'
        else: return string

    def get_flows(self, query=None):
        query_params = self.flow_query_to_query_params(query=query)
        response = self.requests.get(self.urls['flows'], query_params)
        return self.json_raise_for_status(response=response)

    def flow_query_to_query_params(self, query=None):
        self.validate_flow_query(query=query)
        query_params = {}
        for _filter in query.get('filters' or []):
            if _filter['field'] == 'uuid':
                query_params['uuid'] = _filter['field']['value']
        return query_params

    def validate_flow_query(self, query=None):
        valid_fields = ['uuid']
        for _filter in query.get('filters' or []):
            if _filter['field'] not in valid_fields:
                raise Exception(("invalid filter field '{field}'."
                                " Valid fields are: {valid_fields}").format(
                                    field=_filter['field'],
                                    valid_fields=valid_fields))
            if _filter['field'] == 'uuid':
                self.validate_uuid_filter(_filter=_filter)

    def validate_uuid_filter(self, _filter=None):
        valid_operators = ['=']
        if _filter['operator'] not in valid_operators:
            raise Exception(("invalid filter operator '{operator}'."
                             " Valid operators are: {valid_operators}").format(
                                 operator=_filter['operator'],
                                 valid_operators=valid_operators))

    def patch_flow(self, key=None, patches=None):
        flow_url = self.base_url + 'flows/' + key + '/'
        response = self.requests.patch(flow_url, json=patches)
        return self.json_raise_for_status(response=response)

    def create_job(self, job=None):
        response = self.requests.post(self.urls['jobs'], json=job)
        return self.json_raise_for_status(response=response)

    def get_jobs(self, query=None):
        query_params = self.job_query_to_query_params(query=query)
        response = self.requests.get(self.urls['jobs'], query_params)
        return self.json_raise_for_status(response=response)

    def job_query_to_query_params(self, query=None):
        self.validate_job_query(query=query)
        query_params = {}
        for _filter in query.get('filters' or []):
            if _filter['field'] == 'uuid':
                query_params['uuid'] = _filter['field']['value']
        return query_params

    def validate_job_query(self, query=None):
        valid_fields = ['uuid']
        for _filter in query.get('filters' or []):
            if _filter['field'] not in valid_fields:
                raise Exception(("invalid filter field '{field}'."
                                " Valid fields are: {valid_fields}").format(
                                    field=_filter['field'],
                                    valid_fields=valid_fields))
            if _filter['field'] == 'uuid':
                self.validate_uuid_filter(_filter=_filter)

    def patch_job(self, key=None, patches=None):
        job_url = self.base_url + 'jobs/' + key + '/'
        response = self.requests.patch(job_url, json=patches)
        return self.json_raise_for_status(response=response)

    def flush_mc_db(self):
        url = self.urls['flush']
        response = self.requests.get(url)
        return self.json_raise_for_status(response=response)

    def create_queue(self, queue=None):
        response = self.requests.post(self.urls['queues'], json=queue)
        return self.json_raise_for_status(response=response)

    def claim_queue_items(self, queue_key=None, params=None):
        params = params or {}
        url = '{queues_root}{queue_key}/claim_items/'.format(
            queues_root=self.urls['queues'], queue_key=queue_key)
        response = self.requests.post(url, data=json.dumps(params))
        return self.json_raise_for_status(response=response)
