import logging
import requests


class A2G2_Client(object):
    def __init__(self, base_url=None, request_client=None, logger=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests
        self.logger = logger or logging

    def generate_urls(self):
        return {
            'mols': self.base_url + 'mols/',
            'counts': self.base_url + 'counts/',
        }

    def create_mol(self, mol=None):
        response = self.request_client.post(self.urls['mols'], data=mol)
        return self.json_raise_for_status(response=response)

    def json_raise_for_status(self, response=None):
        response.raise_for_status()
        return response.json()

    def get_counts(self):
        response = self.request_client.get(self.urls['counts'])
        return self.json_raise_for_status(response=response)
