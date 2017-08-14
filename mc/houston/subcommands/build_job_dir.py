import json
import tempfile

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def add_arguments(self, parser):
        parser.add_argument(
            '--job_dict',
            help="JSON string for a job dict",
            required=True, type=json.loads
        )

        parser.add_argument(
            '--output_dir',
            help=("Where to create dir. If dir exists, files will be added to"
                  " this dir.")
        )

        parser.add_argument(
            '--build_work_dir_fn',
            help=("Dot_spec to a build_work_dir_fn. If not specified, will try"
                  " to dispatch based on job_type in job_dict."),
            type=self.dot_spec_arg
        )

    def _ensure_parsed_args(self, *args, **kwargs):
        super()._ensure_parsed_args(*args, **kwargs)
        if (
            ('build_work_dir_fn' in self.parsed_args)
            and isinstance(self.parsed_args['build_work_dir_fn'], str)
        ):
            self.parsed_args['build_work_dir_fn'] = self.load_from_dot_spec(
                self.parsed_args['build_work_dir_fn'])

    def _run(self):
        from mc.utils.job_modules.job_dir_builder import JobDirBuilder
        builder = JobDirBuilder()
        output_dir = self._get_output_dir()
        builder.build_job_dir(
            job_dict=self.parsed_args['job_dict'],
            output_dir=output_dir,
            build_work_dir_fn=self.parsed_args.get('build_work_dir_fn')
        )
        return {'job_dir': output_dir}

    def _get_output_dir(self):
        return self.parsed_args.get('output_dir') or tempfile.mkdtemp()
