import requests
from local_job_runner.job_spec_client import JobSpecClient


class MissionControlJobSpecClient(JobSpecClient):
    def __init__(self, base_url):
        self.base_url = base_url

    def fetch_job_specs(self):
        jobs_url = self.base_url + 'jobs/'
        response = requests.get(jobs_url)
        return response.json()

    def claim_job_specs(self, uuids=None):
        claim_jobs_url = self.base_url + 'claim_jobs/'
        response = requests.post(claim_jobs_url, {'uuids': uuids})
        return response.json()

    def update_job_specs(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            job_spec_url = self.base_url + 'jobs/' + _uuid + '/'
            response_for_uuid = requests.patch(job_spec_url, updates_for_uuid)
            results_by_uuid[_uuid] = response_for_uuid.json()
        return results_by_uuid


