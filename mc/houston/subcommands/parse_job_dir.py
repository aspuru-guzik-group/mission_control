from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Parse a job dir"

    def add_arguments(self, parser):
        parser.add_argument(
            '--job_dir',
            help="path to job dir",
            required=True
        )

        parser.add_argument(
            '--parse_work_dir_fn',
            help=("Dot_spec to a parse_work_dir_fn. If not specified, will try"
                  " to dispatch based on job_type from job metadata file"
                  " in job dir."),
            type=self.dot_spec_arg
        )

    def _ensure_parsed_args(self, *args, **kwargs):
        super()._ensure_parsed_args(*args, **kwargs)
        assert self.parsed_args['job_dir'] is not None, 'job_dir is required'
        if (
            ('parse_work_dir_fn' in self.parsed_args)
            and isinstance(self.parsed_args['parse_work_dir_fn'], str)
        ):
            self.parsed_args['parse_work_dir_fn'] = self.load_from_dot_spec(
                self.parsed_args['parse_work_dir_fn'])

    def _run(self):
        from mc.utils.job_modules.job_dir_parser import JobDirParser
        return JobDirParser().parse_job_dir(**self.parsed_args)
