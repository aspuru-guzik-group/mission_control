"""
Base classes for writing commands.
Adapted from django.core.management.base.py
"""
import io
import os
import sys

import mc
from . import color
from .errors import CommandError
from .argument_parser import ArgumentParser


def handle_default_args(args):
    """
    Include any default args that all commands should accept here
    so that ManagementUtility can handle them before searching for
    user commands.
    """
    if args.get('pythonpath'):
        sys.path.insert(0, args['pythonpath'])


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
        if not msg:
            return
        if not isinstance(msg, str):
            msg = str(msg)
        ending = self.ending if ending is None else ending
        if ending and not msg.endswith(ending):
            msg += ending
        style_func = style_func or self.style_func
        self._out.write(style_func(msg))


class BaseCommand(object):
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
       any environment changes requested by args like
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

    def __init__(self, stdout=None, stderr=None, no_color=False, **kwargs):
        self.stdout = OutputWrapper(stdout or sys.stdout)
        self.stderr = OutputWrapper(stderr or sys.stderr)
        if no_color:
            self.style = color.no_style()
        else:
            self.style = color.color_style()
            self.stderr.style_func = self.style.ERROR

    @classmethod
    def run(cls, argv=None):
        argv = argv if argv is not None else sys.argv
        command = cls()
        command.run_from_argv(argv)

    def get_version(self):
        """
        Return the MC version, which should be correct for all built-in
        MC commands. User-supplied commands can override this method to
        return their own version.
        """
        return mc.get_version()

    def create_parser(self, prog_name):
        """
        Create and return the ``ArgumentParser`` which will be used to
        parse the arguments to this command.
        """
        parser_cls = self.get_base_parser_cls()
        parser = parser_cls(self, prog="%s" % (os.path.basename(prog_name)),
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

    def get_base_parser_cls(self): return ArgumentParser

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
        parser = self.create_parser(argv[0])
        parsed_args_ns, unparsed_args = parser.parse_known_args(args=argv[1:])
        parsed_args = vars(parsed_args_ns)
        handle_default_args(parsed_args)
        try:
            self.execute(parsed_args=parsed_args, unparsed_args=unparsed_args)
        except Exception as e:
            if parsed_args.get('traceback') or not isinstance(e, CommandError):
                raise
            self.stderr.write(msg='%s: %s' % (e.__class__.__name__, e))
            sys.exit(1)

    def execute(self, parsed_args=None, unparsed_args=None):
        """Try to execute this command."""
        if parsed_args.get('no_color'):
            self.style = color.no_style()
            self.stderr.style_func = None
        if parsed_args.get('stdout'):
            self.stdout = OutputWrapper(parsed_args['stdout'])
        if parsed_args.get('stderr'):
            self.stderr = OutputWrapper(parsed_args['stderr'],
                                        self.stderr.style_func)
        output = self.handle(parsed_args=parsed_args,
                             unparsed_args=unparsed_args)
        self.stdout.write(msg=output)
        return output

    def handle(self, parsed_args=None, unparsed_args=None):
        """
        The actual logic of the command. Subclasses must implement
        this method.
        """
        raise NotImplementedError('subclasses of BaseCommand must provide a'
                                  'handle() method')
