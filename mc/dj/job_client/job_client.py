import logging
import requests
from jobs.constants import JobStatuses


class MissionControlJobClient(object):

    Statuses = JobStatuses

    def __init__(self, base_url=None, request_client=None, logger=None):
        self.base_url = base_url
        self.urls = self.generate_urls()
        self.request_client = request_client or requests
        self.logger = logger or logging

    def generate_urls(self):
        return {
            'jobs': self.base_url + 'jobs/',
            'claim_jobs': self.base_url + 'claim_jobs/'
        }

    def create_job(self, job_kwargs=None):
        try:
            response = self.request_client.post(self.urls['jobs'], data=job_kwargs)
            if not str(response.status_code).startswith('2'):
                raise Exception("Bad response: %s" % response)
            return response.json()
        except Exception as e:
            msg = ("Client error, request was: {request}, error was:"
                   " '{error}'"
                  ).format(request={'url': self.urls['jobs'],
                                    'data': job_kwargs},
                           error=e)
            self.logger.error(msg)
            raise e

    def fetch_jobs(self, query_params=None):
        if query_params:
            response = self.request_client.get(self.urls['jobs'], query_params)
        else:
            response = self.request_client.get(self.urls['jobs'])
        return response.json()

    def fetch_claimable_jobs(self):
        return self.fetch_jobs(
            query_params={'status': self.Statuses.Pending.name})

    def fetch_job_by_uuid(self, uuid=None):
        fetch_jobs_result = self.fetch_jobs(
            query_params={'uuid': uuid})
        if len(fetch_jobs_result) == 1: return fetch_jobs_result[0]
        else: return None

    def claim_jobs(self, uuids=None):
        response = self.request_client.post(self.urls['claim_jobs'],
                                            {'uuids': uuids})
        return response.json()

    def update_jobs(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            job_url = self.urls['jobs'] + _uuid + '/'
            response = self.request_client.patch(job_url, updates_for_uuid)
            results_by_uuid[_uuid] = response.json()
        return results_by_uuid
