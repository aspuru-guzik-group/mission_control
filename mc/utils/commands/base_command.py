"""
Base classes for writing commands.
Adapted from django.core.management.base.py
"""
import argparse
import io
import os
import sys

import mc
from . import color


class CommandError(Exception):
    """
    Exception class indicating a problem while executing a management
    command.

    If this exception is raised during the execution of a management
    command, it will be caught and turned into a nicely-printed error
    message to the appropriate output stream (i.e., stderr); as a
    result, raising this exception (with a sensible description of the
    error) is the preferred way to indicate that something has gone
    wrong in the execution of a command.
    """
    pass

class CommandParser(argparse.ArgumentParser):
    """
    Customized ArgumentParser class to improve some error messages and prevent
    SystemExit in several occasions, as SystemExit is unacceptable when a
    command is called programmatically.
    """
    def __init__(self, cmd=None, **kwargs):
        self.cmd = cmd
        super().__init__(**kwargs)

    def parse_args(self, args=None, namespace=None):
        # Catch missing argument for a better error message
        if (hasattr(self.cmd, 'missing_args_message') and
                not (args or any(not arg.startswith('-') for arg in args))):
            self.error(self.cmd.missing_args_message)
        return super().parse_args(args, namespace)

    def error(self, message):
        if self.cmd._called_from_command_line:
            super().error(message)
        else:
            raise CommandError("Error: %s" % message)


def handle_default_options(options):
    """
    Include any default options that all commands should accept here
    so that ManagementUtility can handle them before searching for
    user commands.
    """
    if options.pythonpath:
        sys.path.insert(0, options.pythonpath)


class OutputWrapper(io.TextIOBase):
    """
    Wrapper around stdout/stderr
    """
    @property
    def style_func(self):
        return self._style_func

    @style_func.setter
    def style_func(self, style_func):
        if style_func and self.isatty():
            self._style_func = style_func
        else:
            self._style_func = lambda x: x

    def __init__(self, out, style_func=None, ending='\n'):
        self._out = out
        self.style_func = None
        self.ending = ending

    def __getattr__(self, name):
        return getattr(self._out, name)

    def isatty(self):
        return hasattr(self._out, 'isatty') and self._out.isatty()

    def write(self, msg=None, style_func=None, ending=None):
        if not msg: return
        ending = self.ending if ending is None else ending
        if ending and not msg.endswith(ending): msg += ending
        style_func = style_func or self.style_func
        self._out.write(style_func(msg))


class BaseCommand:
    """
    A base class for common command logic.

    Use this class if you want access to all of the mechanisms which
    parse the command-line arguments and work out what code to call in
    response; if you don't need to change any of that behavior,
    consider using one of the subclasses defined in this file.

    If you are interested in overriding/customizing various aspects of
    the command-parsing and -execution behavior, the normal flow works
    as follows:

    1. Something loads the command class and calls its ``run()`` method.

    2. The ``run_from_argv()`` method calls ``create_parser()`` to get
       an ``ArgumentParser`` for the arguments, parses them, performs
       any environment changes requested by options like
       ``pythonpath``, and then calls the ``execute()`` method,
       passing the parsed arguments.

    3. The ``execute()`` method attempts to carry out the command by
       calling the ``handle()`` method with the parsed arguments; any
       output produced by ``handle()`` will be printed to standard
       output.

    4. If ``handle()`` or ``execute()`` raised any exception (e.g.
       ``CommandError``), ``run_from_argv()`` will  instead print an error
       message to ``stderr``.

    Thus, the ``handle()`` method is typically the starting point for
    subclasses; many built-in commands and command types either place
    all of their logic in ``handle()``, or perform some additional
    parsing work in ``handle()`` and then delegate from it to more
    specialized methods as needed.

    Several attributes affect behavior at various steps along the way:

    ``help``
        A short description of the command, which will be printed in
        help messages.
    """
    help = ''

    _called_from_command_line = False
    def __init__(self, stdout=None, stderr=None, no_color=False):
        self.stdout = OutputWrapper(stdout or sys.stdout)
        self.stderr = OutputWrapper(stderr or sys.stderr)
        if no_color: self.style = color.no_style()
        else:
            self.style = color.color_style()
            self.stderr.style_func = self.style.ERROR

    @classmethod
    def run(cls, argv=None):
        if argv is None: argv = sys.argv
        command = cls()
        command.run_from_argv(argv)

    def get_version(self):
        """
        Return the Django version, which should be correct for all built-in
        Django commands. User-supplied commands can override this method to
        return their own version.
        """
        return mc.get_version()

    def create_parser(self, prog_name):
        """
        Create and return the ``ArgumentParser`` which will be used to
        parse the arguments to this command.
        """
        parser = CommandParser(self, prog="%s" % (os.path.basename(prog_name)),
                               description=self.help or None)
        parser.add_argument('--version', action='version',
                            version=self.get_version())
        parser.add_argument(
            '--pythonpath',
            help=('A directory to add to the Python path,'
                  ' e.g. "/home/myproject".')
        )
        parser.add_argument('--traceback', action='store_true',
                            help='Raise on CommandError exceptions')
        parser.add_argument(
            '--no-color', action='store_true', dest='no_color',
            help="Don't colorize the command output.",
        )
        self.add_arguments(parser=parser)
        return parser

    def add_arguments(self, parser=None):
        """
        Entry point for subclassed commands to add custom arguments.
        """
        pass

    def print_help(self, prog_name):
        """
        Print the help message for this command, derived from
        ``self.usage()``.
        """
        parser = self.create_parser(prog_name)
        parser.print_help()

    def run_from_argv(self, argv):
        """
        Set up any environment changes requested (e.g., Python path settings),
         then run this command. If the
        command raises a ``CommandError``, intercept it and print it sensibly
        to stderr. If the ``--traceback`` option is present or the raised
        ``Exception`` is not ``CommandError``, raise it.
        """
        self._called_from_command_line = True
        parser = self.create_parser(argv[0])
        options, unparsed_args = parser.parse_known_args(args=argv[1:])
        cmd_options = vars(options)
        # Move positional args out of options to mimic legacy optparse
        args = cmd_options.pop('args', ())
        handle_default_options(options)
        try: self.execute(unparsed_args, *args, **cmd_options)
        except Exception as e:
            if options.traceback or not isinstance(e, CommandError): raise
            self.stderr.write(msg='%s: %s' % (e.__class__.__name__, e))
            sys.exit(1)

    def execute(self, unparsed_args, *args, **options):
        """Try to execute this command."""
        if options['no_color']:
            self.style = color.no_style()
            self.stderr.style_func = None
        if options.get('stdout'):
            self.stdout = OutputWrapper(options['stdout'])
        if options.get('stderr'):
            self.stderr = OutputWrapper(options['stderr'],
                                        self.stderr.style_func)
        output = self.handle(args=args, kwargs=options,
                             unparsed_args=unparsed_args)
        self.stdout.write(msg=output)
        return output

    def handle(self, args=None, kwargs=None, unparsed_args=None):
        """
        The actual logic of the command. Subclasses must implement
        this method.
        """
        raise NotImplementedError('subclasses of BaseCommand must provide a'
                                  'handle() method')
