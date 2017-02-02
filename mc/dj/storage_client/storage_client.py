import json
import requests
import logging


class MissionControlStorageClient(object):
    def __init__(self, base_url=None, request_client=None, logger=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests
        self.logger = logger or logging

    def generate_urls(self):
        return {
            'post': self.base_url + 'post/',
            'get': self.base_url + 'get/'
        }

    def post_data(self, data=None, storage_params=None):
        storage_params = storage_params or {}
        url = self.urls['post']
        try:
            response = self.request_client.post(url, data={
                'data': data,
                'params': self.serialize_storage_params(storage_params),
            })
            if not str(response.status_code).startswith('2'):
                raise Exception("Bad response: %s" % response)
            result = response.json()
            return result['params']
        except Exception as e:
            msg = ("Client error, request was: {request}, error was:"
                   " '{error}'"
                  ).format(request={'url': url}, error=e)
            self.logger.error(msg)
            raise e

    def serialize_storage_params(self, storage_params=None):
        return json.dumps(storage_params)

    def get_data(self, storage_params=None):
        url = self.urls['get']
        try:
            response = self.request_client.get(url, data={
                'params': self.serialize_storage_params(storage_params),
            })
            if not str(response.status_code).startswith('2'):
                raise Exception("Bad response: %s" % response)
            return response.json()
        except Exception as e:
            msg = ("Client error, request was: {request}, error was:"
                   " '{error}'"
                  ).format(request={'url': url}, error=e)
            self.logger.error(msg)
            raise e
