import os

from .bash import BashSubmissionBuilder


class SlurmSubmissionBuilder(object):
    def __init__(self, *args, submission_spec=None, **kwargs):
        self.bash_builder = BashSubmissionBuilder(*args, **kwargs)

    def build_submission(self, submission_spec=None, output_dir=None, **kwargs):
        self._spec = submission_spec
        self._output_dir = output_dir
        bash_submission_spec = {
            **self._spec,
            'header': "{sbatch_header}\n{header}".format(
                sbatch_header=self.generate_sbatch_header(),
                header=submission_spec.get('header')
            )
        }
        return self.bash_builder.build_submission(
            submission_spec=bash_submission_spec, **kwargs)

    def generate_sbatch_header(self):
        sbatch_params = {
            'error': os.path.join(
                self._output_dir,
                self.bash_builder.std_log_file_names['stderr']
            ),
            'output': os.path.join(
                self._output_dir,
                self.bash_builder.std_log_file_names['stdout']
            )
            **self.spec.get('sbatch_params', {})
        }
        sbatch_lines = [
            "#SBATCH --{key}={value}".format(key=key, value=value)
            for key, value in sbatch_params.items()
        ]
        return "\n".join(sbatch_lines)

