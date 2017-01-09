import requests
from jobs.models import JobStatuses


class MissionControlJobSpecClient(object):

    Statuses = JobStatuses

    def __init__(self, base_url=None, request_client=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests

    def generate_urls(self):
        return {
            'jobs': self.base_url + 'jobs/',
            'claim_jobs': self.base_url + 'claim_jobs/'
        }

    def create_job(self, job=None):
        response = self.request_client.post(self.urls['jobs'], data=job)
        return response.json()

    def fetch_job_specs(self, query_params=None):
        if query_params:
            response = self.request_client.get(self.urls['jobs'], query_params)
        else:
            response = self.request_client.get(self.urls['jobs'])
        return response.json()

    def fetch_claimable_job_specs(self):
        return self.fetch_job_specs(
            query_params={'status': self.Statuses.Pending.name})

    def fetch_job_spec_by_uuid(self, uuid=None):
        fetch_job_specs_result = self.fetch_job_specs(
            query_params={'uuid': uuid})
        if len(fetch_job_specs_result) == 1: return fetch_job_specs_result[0]
        else: return None

    def claim_job_specs(self, uuids=None):
        response = self.request_client.post(self.urls['claim_jobs'],
                                            {'uuids': uuids})
        return response.json()

    def update_job_specs(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            job_spec_url = self.urls['jobs'] + _uuid + '/'
            response = self.request_client.patch(job_spec_url, updates_for_uuid)
            results_by_uuid[_uuid] = response.json()
        return results_by_uuid
