from mc.utils.commands.base_command import BaseCommand


class SubcommandCommand(BaseCommand):
    subcommands = []

    class InvalidSubcommandError(Exception):
        def __init__(self, subcommand=None, **kwargs):
            msg = "subcommand was '%s'" % subcommand
            super().__init__(msg, **kwargs)

    def add_arguments(self, parser=None):
        subparsers = parser.add_subparsers(title='subcommand',
                                           dest='subcommand')
        for subcommand in self.subcommands: subparsers.add_parser(subcommand)

    def handle(self, args=None, kwargs=None, unparsed_args=None):
        subcommand = kwargs.get('subcommand')
        try: subcommand_fn = getattr(self, subcommand)
        except Exception as exc:
            raise self.InvalidSubcommandError(subcommand) from exc
        return subcommand_fn(args=args, kwargs=kwargs,
                             unparsed_args=unparsed_args)
