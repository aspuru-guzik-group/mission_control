class JobEngine(object):
    def __init__(self, *args, command_handlers=None, **kwargs):
        self.command_handlers = command_handlers or {}

    def execute_job(self, *args, job=None, **kwargs):
        command = job['job_spec']['command']
        self.handle_command(*args, command=command, job=job, **kwargs)

    def handle_command(self, *args, command=None, job=None, **kwargs):
        try:
            handler = self.command_handlers[command]
            handler(*args, job=job, **kwargs)
        except KeyError:
            raise Exception("No handler found for command '{command}'".format(
                command=command))
