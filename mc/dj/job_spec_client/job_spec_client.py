import requests
from local_job_runner.job_spec_client import JobSpecClient
from jobs.models import Job


class MissionControlJobSpecClient(JobSpecClient):
    def __init__(self, base_url=None, request_client=None):
        self.base_url = base_url
        self.request_client = request_client or requests

    def fetch_job_specs(self, filters=None):
        jobs_url = self.base_url + 'jobs/'
        if filters:
            response = self.request_client.get(jobs_url, filters)
        else:
            response = self.request_client.get(jobs_url)
        return response.json()

    def fetch_claimable_job_specs(self):
        return self.fetch_job_specs(filters={'status': Job.STATUSES.PENDING})

    def fetch_candidate_job_specs(self):
        jobs_url = self.base_url + 'jobs/'
        response = self.request_client.get(jobs_url)
        return response.json()

    def claim_job_specs(self, uuids=None):
        claim_jobs_url = self.base_url + 'claim_jobs/'
        response = self.request_client.post(claim_jobs_url, {'uuids': uuids})
        return response.json()

    def update_job_specs(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            job_spec_url = self.base_url + 'jobs/' + _uuid + '/'
            response_for_uuid = self.request_client.patch(
                job_spec_url, updates_for_uuid)
            results_by_uuid[_uuid] = response_for_uuid.json()
        return results_by_uuid


