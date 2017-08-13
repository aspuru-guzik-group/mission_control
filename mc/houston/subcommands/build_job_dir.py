import json
import tempfile

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def add_arguments(self, parser):
        parser.add_argument('--job_dict', help="JSON string for a job dict",
                            required=True, type=json.loads)

        parser.add_argument(
            '--output_dir',
            help=(
                "Where to put job files."
                " Dir will be created if it does not exist."
                " If not specified a tmp dir will be used."
            ),
        )

    def _ensure_parsed_args(self, *args, **kwargs):
        super()._ensure_parsed_args(*args, **kwargs)
        if not self.parsed_args.get('output_dir'):
            self.parsed_args['output_dir'] = tempfile.mkdtemp()

    def _run(self):
        from mc.utils.job_modules.job_dir_builder import JobDirBuilder
        builder = JobDirBuilder()
        builder.build_job_dir(
            job_dict=self.parsed_args['job_dict'],
            output_dir=self.parsed_args['output_dir']
        )
        return {'job_dir': self.parsed_args['output_dir']}
