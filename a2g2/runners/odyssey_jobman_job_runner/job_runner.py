import json
import os
import logging
import tempfile


class JobRunner(object):
    def __init__(self, job_client=None, job_submission_factory=None,
                 jobman=None, logging_cfg=None, max_claims_per_tick=None,
                 jobman_source_name=None, **job_runner_kwargs):
        self.job_client = job_client
        self.job_submission_factory = job_submission_factory
        self.jobman = jobman
        self.max_claims_per_tick = max_claims_per_tick or 3
        self.logger = self._generate_logger(logging_cfg=logging_cfg)
        self.jobman_source_name = jobman_source_name or 'mc'

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
            'artifact_spec': self.generate_artifact_spec_for_dir(
                _dir=jobman_job['submission']['dir']),
            'std_log_contents': self.jobman.get_std_log_contents_for_job(
                job=jobman_job),
        }
        parsed_jobman_job['error'] = \
                parsed_jobman_job['std_log_contents'].get('failure')
        parsed_jobman_job['status'] = 'COMPLETED'
        if parsed_jobman_job['jobman_job']['status'] == 'FAILED' or \
           parsed_jobman_job['error']: parsed_jobman_job['status'] = 'FAILED'
        return parsed_jobman_job

    def generate_artifact_spec_for_dir(self, _dir=None):
        artifact_spec = {
            'artifact_type': 'a2g2.artifacts.odyssey',
            'artifact_params': {'path': _dir}
        }
        return artifact_spec

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
                'artifact': parsed_jobman_job.get('artifact_spec'),
                'std_log_contents': self.get_std_log_contents(
                    parsed_jobman_job=parsed_jobman_job)
            }
        }
        return patch

    def get_std_log_contents(self, parsed_jobman_job=None):
        mc_job = parsed_jobman_job['jobman_job']['source_meta']['mc_job']
        contents = parsed_jobman_job['std_log_contents']
        logs_to_expose = mc_job['job_spec'].get('std_logs_to_expose', [])
        if logs_to_expose == 'all':
            contents_to_expose = contents
        else:
            contents_to_expose = {
                log_to_expose: contents.get(log_to_expose)
                for log_to_expose in logs_to_expose
            }
        return contents_to_expose

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
        self.prepare_job_inputs(mc_job=mc_job, submission_dir=submission_dir)
        submission = self.job_submission_factory.build_job_submission(
            job=mc_job, submission_dir=submission_dir)
        return submission

    def prepare_job_inputs(self, mc_job=None, submission_dir=None):
        inputs_dir = os.path.join(submission_dir, 'inputs')
        os.makedirs(inputs_dir, exist_ok=True)
        artifacts = mc_job['job_spec'].get('inputs', {}).get('artifacts', {})
        for artifact_key, artifact in artifacts.items():
            self.prepare_input_artifact(artifact_key=artifact_key,
                                        artifact=artifact,
                                        inputs_dir=inputs_dir)

    def prepare_input_artifact(self, artifact_key=None, artifact=None,
                               inputs_dir=None):
        if artifact['artifact_type'] == 'a2g2.artifacts.odyssey':
            os.symlink(artifact['artifact_params']['path'],
                       os.path.join(inputs_dir, artifact_key))
