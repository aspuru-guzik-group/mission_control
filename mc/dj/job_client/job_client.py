import json
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
            formatted_job_kwargs = self.format_job_kwargs(job_kwargs)
            response = self.request_client.post(self.urls['jobs'],
                                                data=formatted_job_kwargs)
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

    def format_job_kwargs(self, job_kwargs=None):
        formatted_job_kwargs = {**job_kwargs}
        if 'data' in job_kwargs:
            formatted_job_kwargs['data'] = self.serialize_job_data(
                job_kwargs['data'])
        return formatted_job_kwargs

    def serialize_job_data(self, job_data=None):
        return json.dumps(job_data)

    def fetch_jobs(self, query_params=None):
        if query_params:
            response = self.request_client.get(self.urls['jobs'], query_params)
        else:
            response = self.request_client.get(self.urls['jobs'])
        jobs = [self.format_fetched_job(fetched_job)
                for fetched_job in response.json()]
        return jobs

    def format_fetched_job(self, fetched_job=None):
        formatted_job = {**fetched_job}
        if fetched_job.get('data', None):
            formatted_job['data'] = self.deserialize_job_data(
                fetched_job['data'])
        return formatted_job

    def deserialize_job_data(self, serialized_job_data=None):
        return json.loads(serialized_job_data)

    def fetch_claimable_jobs(self):
        return self.fetch_jobs(
            query_params={'status': self.Statuses.PENDING.name})

    def fetch_job_by_uuid(self, uuid=None):
        result = self.fetch_jobs(query_params={'uuid': uuid})
        if len(result) == 1: return result[0]
        else: raise Exception("Could not get job with uuid '%s'" % uuid)

    def claim_jobs(self, uuids=None):
        response = self.request_client.post(self.urls['claim_jobs'],
                                            {'uuids': uuids})
        return response.json()

    def update_jobs(self, updates_by_uuid=None):
        results_by_uuid = {}
        for _uuid, updates_for_uuid in updates_by_uuid.items():
            job_url = self.urls['jobs'] + _uuid + '/'
            formatted_job_kwargs = self.format_job_kwargs(
                job_kwargs=updates_for_uuid)
            response = self.request_client.patch(job_url,
                                                 data=formatted_job_kwargs)
            results_by_uuid[_uuid] = response.json()
        return results_by_uuid
