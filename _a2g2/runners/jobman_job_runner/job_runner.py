import json
import os
import logging
import tempfile


class JobRunner(object):
    def __init__(self, cfg=None, **kwargs):
        self.cfg = cfg
        self.job_client = cfg.get('job_client')
        self.submission_factory = cfg.get('submission_factory')
        self.jobman = cfg.get('jobman')
        self.max_claims_per_tick = cfg.get('max_claims_per_tick') or 3
        self.logger = self._generate_logger(logging_cfg=cfg.get('logging_cfg'))
        self.jobman_source_name = cfg.get('jobman_source_name') or 'mc'
        self.artifact_processor = cfg.get('artifact_processor')
        self.get_cfg_for_mc_job = cfg.get('get_cfg_for_mc_job')  or \
                self._default_get_cfg_for_mc_job

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
    
    def _default_get_cfg_for_mc_job(*args, **kwargs): return {}

    def tick(self):
        self.logger.info("tick")
        self.process_executed_jobs()
        self.fill_jobman_queue()

    def process_executed_jobs(self):
        executed_jobman_jobs = self.get_executed_jobman_jobs()
        parsed_jobman_jobs = self.parse_jobman_jobs(
            jobman_jobs=executed_jobman_jobs)
        keyed_patches = self.parsed_jobman_jobs_to_keyed_patches(
            parsed_jobman_jobs=parsed_jobman_jobs)
        self.patch_jobs(keyed_patches=keyed_patches)
        self.finalize_jobman_jobs(jobman_jobs=executed_jobman_jobs)

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
            'artifact': self.artifact_processor.dir_to_artifact(
                _dir=jobman_job['submission']['dir']),
            'std_log_contents': self.get_std_log_contents_for_jobman_job(
                jobman_job=jobman_job),
        }
        parsed_jobman_job['error'] = \
                parsed_jobman_job['std_log_contents'].get('failure')
        parsed_jobman_job['status'] = 'COMPLETED'
        if parsed_jobman_job['jobman_job']['status'] == 'FAILED' or \
           parsed_jobman_job['error']: parsed_jobman_job['status'] = 'FAILED'
        return parsed_jobman_job

    def get_std_log_contents_for_jobman_job(self, jobman_job=None):
        submission = jobman_job['submission']
        mc_job = jobman_job['source_meta']['mc_job']
        logs_to_expose = mc_job['job_spec'].get('std_logs_to_expose', [])
        if logs_to_expose == 'all':
            logs_to_expose = submission['std_log_files'].keys()
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
        rel_log_path = submission['std_log_files'][log]
        abs_log_path = os.path.join(submission['dir'], rel_log_path)
        if os.path.exists(abs_log_path):
            with open(abs_log_path) as f: return f.read()

    def parsed_jobman_jobs_to_keyed_patches(self, parsed_jobman_jobs=None):
        keyed_patches = {}
        for parsed_jobman_job in parsed_jobman_jobs:
            patch = self.parsed_jobman_job_to_patch(
                parsed_jobman_job=parsed_jobman_job)
            mc_job = parsed_jobman_job['jobman_job']['source_meta']['mc_job']
            keyed_patches[mc_job['uuid']] = patch
        return keyed_patches

    def parsed_jobman_job_to_patch(self, parsed_jobman_job=None):
        patch = {
            'status': parsed_jobman_job['status'],
            'data': {
                'artifact': parsed_jobman_job.get('artifact'),
                'std_log_contents': parsed_jobman_job.get('std_log_contents')
            }
        }
        return patch

    def deserialize_mc_job(self, mc_job=None):
        deserialized_mc_job = {}
        for k, v in mc_job.items():
            if k == 'data': v = json.loads(v or '{}')
            deserialized_mc_job[k] = v
        return deserialized_mc_job

    def patch_jobs(self, keyed_patches=None):
        keyed_serialized_patches = {
            key: self.serialize_job_patch(patch=patch)
            for key, patch in keyed_patches.items()
        }
        self.job_client.patch_jobs(keyed_patches=keyed_serialized_patches)

    def serialize_job_patch(self, patch=None):
        serialized_patch = {}
        for k, v in (patch or {}).items():
            if k == 'data': v = json.dumps(v or '{}')
            serialized_patch[k] = v
        return serialized_patch

    def finalize_jobman_jobs(self, jobman_jobs=None):
        finalized_jobman_jobs = [{**jobman_job, 'status': 'COMPLETED'}
                                 for jobman_job in jobman_jobs]
        self.jobman.save_jobs(jobs=finalized_jobman_jobs)

    def fill_jobman_queue(self):
        for mc_job in self.job_client.claim_jobs(params={
            'limit': self.get_claim_limit()
        }):
            self.jobman.submit_job(
                submission=self.build_submission(mc_job=mc_job),
                source=self.jobman_source_name,
                source_meta={'mc_job': mc_job}
            )

    def get_claim_limit(self):
        return min(self.max_claims_per_tick, self.jobman.num_free_slots)

    def build_submission(self, mc_job=None):
        submission_dir = tempfile.mkdtemp(prefix='sf.')
        cfg_for_mc_job = self.get_cfg_for_mc_job(mc_job=mc_job)
        self.prepare_job_inputs(mc_job=mc_job, cfg=cfg_for_mc_job,
                                submission_dir=submission_dir)
        submission = self.submission_factory.build_submission(
            job=mc_job, cfg=cfg_for_mc_job, output_dir=submission_dir)
        return submission

    def prepare_job_inputs(self, mc_job=None, cfg=None, submission_dir=None):
        inputs_dir = os.path.join(submission_dir, 'inputs')
        os.makedirs(inputs_dir, exist_ok=True)
        artifacts = mc_job['job_spec'].get('inputs', {}).get('artifacts', {})
        for artifact_key, artifact in artifacts.items():
            dest = os.path.join(inputs_dir, artifact_key)
            self.artifact_processor.artifact_to_dir(
                artifact=artifact, dest=dest)
