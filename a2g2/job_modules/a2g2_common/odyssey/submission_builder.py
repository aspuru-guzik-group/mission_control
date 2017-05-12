import json
import os
import textwrap

from mc.a2g2.job_dir_builders.odyssey import OdysseyJobDirBuilder


class SubmissionBuilder(object):
    def __init__(self, job=None, cfg=None, submission_dir=None):
        self.job = job
        self.cfg = cfg
        self.submission_dir = submission_dir
        self.submission_meta_file_name = 'submission.json'

    def build_submission(self):
        self.ensure_dir(dir=self.submission_dir)
        submission_meta = OdysseyJobDirBuilder.build_dir(
            dir_spec={
                'entrypoint_body': self.generate_entrypoint_body(),
            },
            submission_dir=self.submission_dir,
            submission_meta_file_name=self.submission_meta_file_name
        )
        return submission_meta

    def ensure_dir(self, dir):
        os.makedirs(self.submission_dir, exist_ok=True)

    def generate_entrypoint_body(self):
        job_engine_cfg = self.cfg.get('job_engine', {})
        entrypoint_body = textwrap.dedent(
            """
            {job_engine_preamble}
            {job_engine_exe} \\
                run_job_submission {job_cli_params} {submission_cli_params}
            """
        ).strip().format(
            job_engine_preamble=job_engine_cfg.get('entrypoint_preamble', ''),
            job_engine_exe=job_engine_cfg['job_engine_exe'],
            job_cli_params=self.params_to_cli_args(
                params=self.write_json_params()),
            submission_cli_params=self.params_to_cli_args(
                params={'submission': self.submission_meta_file_name})
        )
        return entrypoint_body

    def write_json_params(self):
        json_params = {}
        for param in ['job', 'cfg']:
            json_filename = param + '.json'
            json_path = os.path.join(self.submission_dir, json_filename)
            with open(json_path, 'w') as f: json.dump(getattr(self, param), f)
            json_params[param] = json_filename
        return json_params

    def params_to_cli_args(self, params=None):
        return ' '.join([
            '--{param}="{value}"'.format(param=param, value=value)
            for param, value in params.items()
        ])

def build_job_submission(*args, job=None, cfg=None, submission_dir=None, **kwargs):
    builder = SubmissionBuilder(job=job, cfg=cfg, submission_dir=submission_dir)
    return builder.build_submission()
