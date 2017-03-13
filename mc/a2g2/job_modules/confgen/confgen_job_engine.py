from ..utils.job_engine import JobEngine
from . import confgen_generator
from . import confgen_parser


class ConformerJobEngine(JobEngine):
    def __init__(self, *args, **kwargs):
        command_handlers = {
            method.__name__: method for method in [
                self.generate_conformers,
                self.parse_completed_confgen_dir,
            ]
        }
        kwargs = {'command_handlers': command_handlers, **kwargs}
        super().__init__(self, *args, **kwargs)

    def generate_conformers(self, *args, job=None, output_dir=None, **kwargs):
        confgen_generator.generate_conformers(**{
            **job['job_spec']['kwargs'],
            'output_dir': output_dir
        })

    def parse_completed_confgen_dir(self, *args, job=None, output_dir=None,
                                    **kwargs):
        confgen_parser.parse_completed_confgen_dir(
            completed_confgen_dir=job['data']['input']['dir_to_parse'],
            output_dir=output_dir
        )

def generate_job_engine():
    return ConformerJobEngine()
