from pathlib import Path
import textwrap

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = 'Create a MissionControl project dir.'

    JOB_DIR_ROOT_REL_PATH = 'job_dirs'
    CFG_REL_PATH = 'houston_cfg.py'
    DB_REL_PATH = 'mc.db.sqlite'

    def add_arguments(self, parser):
        parser.add_argument(
            '--output_path',
            help="path at which to create dir",
            required=True
        )

    def _run(self):
        output_path = Path(self.parsed_args['output_path'])
        output_path.mkdir(parents=True)
        (output_path / self.JOB_DIR_ROOT_REL_PATH).mkdir()
        cfg_path = (output_path / self.CFG_REL_PATH)
        with cfg_path.open('w') as f:
            f.write(self._generate_cfg_content())

    def _generate_cfg_content(self):
        return textwrap.dedent(
            '''
            from pathlib import Path

            dir_path = Path(__file__).parent.absolute()

            MC_DB_URI = 'sqlite:///%s' % (dir_path / '{DB_REL_PATH}')

            JOB_DIRS_ROOT = str(dir_path / '{JOB_DIR_ROOT_REL_PATH}')
            '''
        ).lstrip().format(
            DB_REL_PATH=self.DB_REL_PATH,
            JOB_DIR_ROOT_REL_PATH=self.JOB_DIR_ROOT_REL_PATH
        )
