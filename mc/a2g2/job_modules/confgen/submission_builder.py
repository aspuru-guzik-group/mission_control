import json
import os
import textwrap

from a2g2.job_dir_builders.odyssey import OdysseyJobDirBuilder


class SubmissionBuilder(object):
    def __init__(self, job=None, cfg=None, submission_dir=None):
        self.job = job
        self.cfg = cfg
        self.submission_dir = submission_dir

    def build_submission(self):
        self.ensure_dir(dir=self.submission_dir)
        submission_meta = OdysseyJobDirBuilder.build_dir(
            dir_spec={
                'entrypoint_body': self.generate_entrypoint_body(),
            },
            output_dir=self.submission_dir
        )
        return submission_meta

    def ensure_dir(self, dir):
        os.makedirs(self.submission_dir, exist_ok=True)

    def generate_entrypoint_body(self):
        job_engine_cfg = self.cfg.get('job_engine', {})
        entrypoint_body = textwrap.dedent(
            """
            {job_engine_preamble}
            python -m {job_engine_module} {job_engine_command} {params}
            """
        ).strip().format(
            job_engine_preamble=job_engine_cfg.get('entrypoint_preamble', ''),
            job_engine_module=job_engine_cfg['engine_module'],
            job_engine_command='run_job_submission',
            params=self.params_to_cli_args(params=self.write_json_params())
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
