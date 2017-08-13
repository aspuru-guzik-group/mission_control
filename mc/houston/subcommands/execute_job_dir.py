from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Run a job_dir. Normally only used for faking execution in tests."

    def add_arguments(self, parser):
        defaults = self._get_defaults()
        parser.add_argument('--job_dir', required=True)
        parser.add_argument('--fake', action='store_true',
                            default=defaults['fake'])

    def _get_defaults(self):
        return {
            'fake': False,
        }

    def _run(self):
        from mc.utils.job_modules.job_dir_runner import JobDirRunner
        runner = JobDirRunner()
        runner.run_job_dir(
            job_dir=self.parsed_args['job_dir'],
            fake=self.parsed_args['fake']
        )
