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

    PROCESSED_TAG = 'PROCESSED'

    def __init__(self, job_record_client=None, build_jobdir_fn=None,
                 jobman=None, max_claims_per_tick=None, logger=None,
                 logging_cfg=None, jobman_source_name=None,
                 artifact_handler=None, jobdirs_dir=None, **kwargs):
        self.logger = logger or self._generate_logger(logging_cfg=logging_cfg)
        self.job_record_client = job_record_client
        self.build_jobdir_fn = build_jobdir_fn
        self.jobman = jobman
        self.max_claims_per_tick = max_claims_per_tick or 3
        self.jobman_source_name = jobman_source_name or 'mc'
        self.artifact_handler = artifact_handler
        self.jobdirs_dir = jobdirs_dir

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
        finished_stats = self.process_finished_jobs()
        claimed_stats = self.fill_jobman_queue()
        tick_stats = {**finished_stats, **claimed_stats}
        return tick_stats

    def process_finished_jobs(self):
        unprocessed_finished_jobman_jobs = \
                self.get_unprocessed_finished_jobman_jobs()
        parsed_jobman_jobs = self.parse_jobman_jobs(
            jobman_jobs=unprocessed_finished_jobman_jobs)
        keyed_patches = self.parsed_jobman_jobs_to_keyed_patches(
            parsed_jobman_jobs=parsed_jobman_jobs)
        self.patch_job_records(keyed_patches=keyed_patches)
        self.finalize_jobman_jobs(jobman_jobs=unprocessed_finished_jobman_jobs)
        finished_stats = {'finished': len(unprocessed_finished_jobman_jobs)}
        return finished_stats

    def get_unprocessed_finished_jobman_jobs(self):
        return self.jobman.get_jobs(query={
            'filters': [
                {'field': 'status', 'op': 'IN', 'arg': ['COMPLETED', 'FAILED']},
                {'field': 'source', 'op': '=', 'arg': self.jobman_source_name},
                {'field': 'source_tag', 'op': '! =', 'arg': self.PROCESSED_TAG},
            ]
        })

    def parse_jobman_jobs(self, jobman_jobs=None):
        return [self.parse_jobman_job(jobman_job=jobman_job)
                for jobman_job in jobman_jobs]

    def parse_jobman_job(self, jobman_job=None):
        parsed_jobman_job = {
            'jobman_job': jobman_job,
            'mc_job': jobman_job['source_meta']['mc_job'],
            'artifact': self.artifact_handler.dir_to_artifact(
                dir_=jobman_job['job_spec']['dir']),
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
        job_spec = jobman_job['job_spec']
        mc_job = jobman_job['source_meta']['mc_job']
        logs_to_expose = mc_job.get('cfg', {}).get('std_logs_to_expose') or []
        if logs_to_expose == 'all':
            logs_to_expose = job_spec.get('std_log_file_names', {}).keys()
        std_log_contents = self.read_jobdir_logs(job_spec=job_spec,
                                                 logs=logs_to_expose)
        return std_log_contents

    def read_jobdir_logs(self, job_spec=None, logs=None):
        logs = logs or []
        return {log: self.read_jobdir_log(job_spec=job_spec, log=log)
                for log in logs}

    def read_jobdir_log(self, job_spec=None, log=None):
        rel_log_path = job_spec['std_log_file_names'][log]
        abs_log_path = os.path.join(job_spec['dir'], rel_log_path)
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
        marked_jobman_jobs = [
            {**jobman_job, 'source_tag': self.PROCESSED_TAG, 'purgeable': 1}
            for jobman_job in jobman_jobs
        ]
        self.jobman.save_jobs(jobs=marked_jobman_jobs)

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
        return min(self.max_claims_per_tick, self.jobman.get_num_free_slots())

    def submit_mc_job(self, mc_job=None):
        self.jobman.submit_jobdir(
            job_spec=self.build_jobdir(mc_job=mc_job),
            source=self.jobman_source_name,
            source_meta={'mc_job': mc_job}
        )

    def build_jobdir(self, mc_job=None):
        jobdir = tempfile.mkdtemp(prefix=self.get_jobdir_prefix(mc_job=mc_job),
                                  dir=self.jobdirs_dir)
        self.prepare_job_inputs(mc_job=mc_job, jobdir=jobdir)
        job_spec = self.build_jobdir_fn(job=mc_job, output_dir=jobdir)
        return job_spec

    def get_jobdir_prefix(self, mc_job=None):
        return 'sf.{job_type}.{now}.'.format(
            job_type=mc_job['job_type'],
            now=datetime.datetime.now().isoformat()
        )

    def prepare_job_inputs(self, mc_job=None, jobdir=None):
        inputs_dir = os.path.join(jobdir, 'inputs')
        os.makedirs(inputs_dir, exist_ok=True)
        artifacts = mc_job.get('job_inputs', {}).get('artifacts') or  {}
        for artifact_key, artifact in artifacts.items():
            dest = os.path.join(inputs_dir, artifact_key)
            self.artifact_handler.artifact_to_dir(artifact=artifact, dest=dest)

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
