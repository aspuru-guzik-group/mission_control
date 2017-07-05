import datetime
import json
import os
import logging
import tempfile
import traceback


class JobRunner(object):
    class SubmissionError(Exception):
        def __init__(self, *args, cause=None, mc_job=None, **kwargs):
            error = "cause: {cause}; job_record: {json_mc_job}".format(
                cause=cause, json_mc_job=json.dumps(mc_job, indent=2))
            super().__init__(self, error, *args, **kwargs)

    def __init__(self, job_record_client=None, submission_factory=None,
                 jobman=None, max_claims_per_tick=None, logger=None,
                 logging_cfg=None, jobman_source_name=None,
                 artifact_processor=None, submissions_dir=None, **kwargs):
        self.logger = logger or self._generate_logger(logging_cfg=logging_cfg)
        self.job_record_client = job_record_client
        self.submission_factory = submission_factory
        self.jobman = jobman
        self.max_claims_per_tick = max_claims_per_tick or 3
        self.jobman_source_name = jobman_source_name or 'mc'
        self.artifact_processor = artifact_processor
        self.submissions_dir = submissions_dir

        self.tick_counter = 0

    def _generate_logger(self, logging_cfg=None):
        logging_cfg = logging_cfg or {}
        logger = logging_cfg.get('logger') \
                or logging.getLogger(logging_cfg.get('logger_name'))
        log_file = logging_cfg.get('log_file')
        if log_file: logger.addHandler(
            logging.FileHandler(os.path.expanduser(log_file)))
        level = logging_cfg.get('level')
        if level:
            logger.setLevel(getattr(logging, level))
        fmt = logging_cfg.get('fmt') or \
                '| %(asctime)s | %(name)s | %(levelname)s | %(message)s\n'
        formatter = logging.Formatter(fmt)
        for handler in logger.handlers: handler.setFormatter(formatter)
        return logger

    def tick(self):
        self.tick_counter += 1
        self.logger.debug('{self}, tick #{tick_counter}'.format(
            self=self, tick_counter=self.tick_counter))
        executed_stats = self.process_executed_jobs()
        claimed_stats = self.fill_jobman_queue()
        tick_stats = {**executed_stats, **claimed_stats}
        return tick_stats

    def process_executed_jobs(self):
        executed_jobman_jobs = self.get_executed_jobman_jobs()
        parsed_jobman_jobs = self.parse_jobman_jobs(
            jobman_jobs=executed_jobman_jobs)
        keyed_patches = self.parsed_jobman_jobs_to_keyed_patches(
            parsed_jobman_jobs=parsed_jobman_jobs)
        self.patch_job_records(keyed_patches=keyed_patches)
        self.finalize_jobman_jobs(jobman_jobs=executed_jobman_jobs)
        executed_stats = {'executed': len(executed_jobman_jobs)}
        return executed_stats

    def get_executed_jobman_jobs(self):
        self.jobman.update_job_engine_states(query={
            'filters': [
                {'field': 'source', 'operator': '=',
                 'value': self.jobman_source_name},
            ]
        })
        return self.jobman.get_jobs(query={
            'filters': [
                {'field': 'status', 'operator': 'IN',
                 'value': ['EXECUTED', 'FAILED']},
                {'field': 'source', 'operator': '=',
                 'value': self.jobman_source_name},
            ]
        })

    def parse_jobman_jobs(self, jobman_jobs=None):
        return [self.parse_jobman_job(jobman_job=jobman_job)
                for jobman_job in jobman_jobs]

    def parse_jobman_job(self, jobman_job=None):
        parsed_jobman_job = {
            'jobman_job': jobman_job,
            'mc_job': jobman_job['source_meta']['mc_job'],
            'artifact': self.artifact_processor.dir_to_artifact(
                dir_=jobman_job['submission']['dir']),
            'std_logs': self.get_std_log_contents_for_jobman_job(
                jobman_job=jobman_job),
        }
        parsed_jobman_job['error'] = \
                parsed_jobman_job['std_logs'].get('failure')
        parsed_jobman_job['status'] = 'COMPLETED'
        if parsed_jobman_job['jobman_job']['status'] == 'FAILED' or \
           parsed_jobman_job['error']: parsed_jobman_job['status'] = 'FAILED'
        return parsed_jobman_job

    def get_std_log_contents_for_jobman_job(self, jobman_job=None):
        submission = jobman_job['submission']
        mc_job = jobman_job['source_meta']['mc_job']
        logs_to_expose = mc_job['job_spec'].get('std_logs_to_expose', [])
        if logs_to_expose == 'all':
            logs_to_expose = submission.get('std_log_file_names', {}).keys()
        std_log_contents = self.read_submission_logs(submission=submission,
                                                     logs=logs_to_expose)
        return std_log_contents

    def read_submission_logs(self, submission=None, logs=None):
        logs = logs or []
        return {
            log: self.read_submission_log(submission=submission, log=log)
            for log in logs
        }

    def read_submission_log(self, submission=None, log=None):
        rel_log_path = submission['std_log_file_names'][log]
        abs_log_path = os.path.join(submission['dir'], rel_log_path)
        if os.path.exists(abs_log_path):
            with open(abs_log_path) as f: return f.read()

    def parsed_jobman_jobs_to_keyed_patches(self, parsed_jobman_jobs=None):
        keyed_patches = {}
        for parsed_jobman_job in parsed_jobman_jobs:
            patch = self.parsed_jobman_job_to_patch(
                parsed_jobman_job=parsed_jobman_job)
            keyed_patches[parsed_jobman_job['mc_job']['key']] = patch
        return keyed_patches

    def parsed_jobman_job_to_patch(self, parsed_jobman_job=None):
        patch = {
            'status': parsed_jobman_job['status'],
            'data': {
                **(parsed_jobman_job['mc_job'].get('data', {})),
                'artifact': parsed_jobman_job.get('artifact'),
                'std_logs': parsed_jobman_job.get('std_logs')
            }
        }
        return patch

    def patch_job_records(self, keyed_patches=None):
        if not keyed_patches: return {}
        self.job_record_client.patch_job_records(keyed_patches=keyed_patches)

    def finalize_jobman_jobs(self, jobman_jobs=None):
        if not jobman_jobs: return
        finalized_jobman_jobs = [{**jobman_job, 'status': 'COMPLETED'}
                                 for jobman_job in jobman_jobs]
        self.jobman.save_jobs(jobs=finalized_jobman_jobs)

    def fill_jobman_queue(self):
        limit = self.get_claim_limit()
        claimed_job_records = self.job_record_client.claim_job_records(
            params={'limit': limit})
        job_records_by_submission_outcome = {'submitted': [],'failed': []}
        for mc_job in claimed_job_records:
            try:
                self.submit_mc_job(mc_job=mc_job)
                job_records_by_submission_outcome['submitted'].append(mc_job)
            except Exception as exc:
                self.logger.exception("SubmissionError")
                error = self.SubmissionError(cause=traceback.format_exc(),
                                             mc_job=mc_job)
                mc_job['data']['error'] = str(error)
                job_records_by_submission_outcome['failed'].append(mc_job)
        self.patch_job_records_per_submission_outcome(
            job_records_by_submission_outcome=job_records_by_submission_outcome)
        claimed_stats = {
            'claimed': len(claimed_job_records),
            **{
                outcome: len(job_records)
                for outcome, job_records in (
                    job_records_by_submission_outcome.items())
            }
        }
        return claimed_stats

    def get_claim_limit(self):
        return min(self.max_claims_per_tick, self.jobman.num_free_slots)

    def submit_mc_job(self, mc_job=None):
        self.jobman.submit_job(
            submission=self.build_job_submission(mc_job=mc_job),
            source=self.jobman_source_name,
            source_meta={'mc_job': mc_job}
        )

    def build_job_submission(self, mc_job=None):
        submission_dir = tempfile.mkdtemp(
            prefix=self.get_submission_dir_prefix(mc_job=mc_job),
            dir=self.submissions_dir
        )
        self.prepare_job_inputs(mc_job=mc_job, submission_dir=submission_dir)
        submission = self.submission_factory.build_job_submission(
            job=mc_job, output_dir=submission_dir)
        return submission

    def get_submission_dir_prefix(self, mc_job=None):
        return 'sf.{job_type}.{now}.'.format(
            job_type=mc_job["job_spec"]["job_type"],
            now=datetime.datetime.now().isoformat()
        )

    def prepare_job_inputs(self, mc_job=None, submission_dir=None):
        inputs_dir = os.path.join(submission_dir, 'inputs')
        os.makedirs(inputs_dir, exist_ok=True)
        artifacts = mc_job['job_spec'].get('inputs', {}).get('artifacts', {})
        for artifact_key, artifact in artifacts.items():
            dest = os.path.join(inputs_dir, artifact_key)
            self.artifact_processor.artifact_to_dir(
                artifact=artifact, dest=dest)

    def patch_job_records_per_submission_outcome(
        self, job_records_by_submission_outcome=None):
        keyed_patches = {
            **{
                job_record['key']: {'status': 'RUNNING'}
                for job_record in job_records_by_submission_outcome['submitted']
            },
            **{
                job_record['key']: {'status': 'FAILED',
                                    'data': job_record['data']}
                for job_record in job_records_by_submission_outcome['failed']
            },
        }
        self.patch_job_records(keyed_patches=keyed_patches)
