import json
import sys

from ...a2g2_common.base_submission_runner import BaseSubmissionRunner
from .. import utils as a2g2_db_utils


class SubmissionRunner(BaseSubmissionRunner):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run_submission(self):
        query_params = self.job['job_spec']['job_params']['query_params']
        query_results = self.execute_query(query_params=query_params)
        self.write_query_results(query_results=query_results)

    def execute_query(self, query_params=None):
        a2g2_client = a2g2_db_utils.generate_a2g2_client(cfg=self.cfg)
        return a2g2_client.execute_query(**query_params)

    def write_query_results(self, query_results=None):
        sys.stdout.write(json.dumps(query_results))

def run_job_submission(*args, job=None, cfg=None, submission=None, **kwargs):
    runner = SubmissionRunner(job=job, cfg=cfg, submission=submission)
    return runner.run_submission()
