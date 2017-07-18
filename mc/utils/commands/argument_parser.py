import argparse

from .errors import CommandError


class ArgumentParser(argparse.ArgumentParser):
    """
    Customized ArgumentParser class to improve some error messages and prevent
    SystemExit in several occasions, as SystemExit is unacceptable when a
    command is called programmatically.
    """
    def __init__(self, cmd=None, formatter_class=None, **kwargs):
        self.cmd = cmd
        formatter_class = formatter_class or \
                argparse.ArgumentDefaultsHelpFormatter
        super().__init__(formatter_class=formatter_class, **kwargs)

    def parse_args(self, args=None, namespace=None):
        # Catch missing argument for a better error message
        if (hasattr(self.cmd, 'missing_args_message') and
                not (args or any(not arg.startswith('-') for arg in args))):
            self.error(self.cmd.missing_args_message)
        return super().parse_args(args, namespace)

    def error(self, message):
        if getattr(self.cmd, '_called_from_command_line', None):
            super().error(message)
        else: raise CommandError("Error: %s" % message)
