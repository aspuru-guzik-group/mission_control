import logging


class JobmanExecutionClient(object):
    def __init__(self, jobman=None, logger=None):
        self.jobman = jobman
        self.logger = logger or logging

    def start_execution(self, submission=None):
        job = self.jobman.submit_job(submission=submission)
        execution_meta = {'job_key': job['job_key']}
        return execution_meta

    def get_execution_state(self, execution_meta=None):
        job = self.jobman.get_jobs(job_keys=[execution_meta['job_key']])
        execution_state = {
            'run_status': job['status'],
            'jobman_job': job
        }
        return execution_state
