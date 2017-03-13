from mc.a2g2.utils import base_command as a2g2_base_command


class BaseCommand(a2g2_base_command.BaseCommand):
    help = 'base_command'

    def add_arguments(self, parser=None):
        super().add_arguments(parser=parser)
        parser.add_argument('--job', type=self.json_file_type)
        parser.add_argument('--ctx_dir', type=str)
        parser.add_argument('--output_dir', type=str)

    def handle(self, *args, **kwargs):
        self.execute_job(*args, **kwargs)

    def execute_job(self, *args, job=None, cfg=None, ctx_dir=None, 
                    output_dir=None, **kwargs): pass
