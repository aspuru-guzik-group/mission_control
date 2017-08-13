import json
import tempfile

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def add_arguments(self, parser):
        parser.add_argument('--job_dict', help="JSON string for a job dict",
                            required=True, type=json.loads)

        parser.add_argument('--output_dir',
                            help=("Where to create dir. If dir exists, files"
                                  " will be added to this dir."))

    def _run(self):
        from a2g2_v2.utils.job_modules.job_dir_builder import JobDirBuilder
        builder = JobDirBuilder()
        output_dir = self._get_output_dir()
        builder.build_job_dir(job_dict=self.parsed_args['job_dict'],
                              output_dir=output_dir)
        return {'job_dir': output_dir}

    def _get_output_dir(self):
        return self.parsed_args.get('output_dir') or tempfile.mkdtemp()
