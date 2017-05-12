from ..flow.submission_runner import SubmissionRunner as FlowSubmissionRunner


class SubmissionRunner(FlowSubmissionRunner):
    def get_flow_spec(self):
        job_params = self.job['job_spec']['job_params']
        flow_spec = {
            'data': job_params.get('flow_data'),
            'node_specs': [
                {'node': {'node_tasks': job_params['tasks']},
                 'precursor_keys': ['ROOT']}
            ]
        }
        return flow_spec

def run_job_submission(*args, job=None, cfg=None, submission=None, **kwargs):
    runner = SubmissionRunner(job=job, cfg=cfg, submission=submission)
    return runner.run_submission()
