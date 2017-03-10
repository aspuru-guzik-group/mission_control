from . import conformer_generator


class ConfgenJobEngine(object):
    def execute_job(self, *args, job=None, **kwargs):
        command = job['job_spec']['command']
        self.dispatch_to_command(*args, command=command, job=job, **kwargs)

    def dispatch_to_command(self, *args, command=None, job=None, **kwargs):
        try:
            getattr(self, command)(*args, job=job, **kwargs)
        except AttributeError:
            raise Exception("No handler found for command '{command}'".format(
                command=command))

    def generate_conformers(self, *args, job=None, output_dir=None, **kwargs):
        conformer_generator.generate_conformers(**{
            **job['job_spec']['kwargs'],
            'output_dir': output_dir
        })

default_job_engine = ConfgenJobEngine()
