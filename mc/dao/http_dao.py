import requests as _requests
import logging

from .base_dao import BaseDao

class HttpDao(BaseDao):
    def __init__(self, base_url=None, urls=None, requests=None, logger=None,
                 max_error_len=1024):
        self.logger = logger or logging
        self.base_url = base_url
        self.urls = urls or self.generate_urls()
        self.requests = requests or _requests
        self.max_error_len = max_error_len

    def generate_urls(self):
        return {
            'Flow': self.base_url + 'flows/',
            'flush': self.base_url + 'flush/',
            'Job': self.base_url + 'jobs/',
            'Queue': self.base_url + 'queues/',
        }

    def create_item(self, item_type=None, kwargs=None):
        response = self.requests.post(self.urls[item_type], json=kwargs)
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

    def get_items(self, item_type=None, query=None):
        query_params = self.query_to_query_params(query=query)
        response = self.requests.get(self.urls[item_type], query_params)
        return self.json_raise_for_status(response=response)

    def query_to_query_params(self, query=None):
        self.validate_query(query=query)
        query_params = {}
        for _filter in query.get('filters' or []):
            if _filter['field'] == 'uuid':
                query_params['uuid'] = _filter['field']['value']
        return query_params

    def patch_item(self, item_type=None, key=None, patches=None):
        url = self.urls[item_type] + key + '/'
        response = self.requests.patch(url, json=patches)
        return self.json_raise_for_status(response=response)

    def flush_mc_db(self):
        url = self.urls['flush']
        response = self.requests.get(url)
        return self.json_raise_for_status(response=response)

    def claim_queue_items(self, queue_key=None):
        url = '{queues_root}{queue_key}/claim_items/'.format(
            queues_root=self.urls['Queue'], queue_key=queue_key)
        response = self.requests.post(url)
        return self.json_raise_for_status(response=response)
