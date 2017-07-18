from .base_command import BaseCommand


class SubcommandCommand(BaseCommand):
    subcommands = []

    class InvalidSubcommandError(Exception):
        def __init__(self, subcommand=None, **kwargs):
            msg = "subcommand was '%s'" % subcommand
            super().__init__(msg, **kwargs)

    class SubcommandExecutionError(Exception):
        def __init__(self, subcommand=None, **kwargs):
            msg = "subcommand was '%s'" % subcommand
            super().__init__(msg, **kwargs)

    def add_arguments(self, parser=None):
        subparsers = parser.add_subparsers(title='subcommand',
                                           dest='subcommand')
        subparsers.required = True
        for subcommand in self.subcommands: subparsers.add_parser(subcommand)

    def handle(self, parsed_args=None, unparsed_args=None):
        self._run_subcommand(subcommand=parsed_args['subcommand'],
                             parsed_args=parsed_args,
                             unparsed_args=unparsed_args)

    def _run_subcommand(self, subcommand=None, parsed_args=None,
                        unparsed_args=None):
        try: subcommand_fn = self._get_subcommand_fn(subcommand=subcommand)
        except Exception as exc:
            raise self.InvalidSubcommandError(subcommand) from exc
        try:
            return self._call_subcommand_fn(subcommand_fn=subcommand_fn,
                                            parsed_args=parsed_args,
                                            unparsed_args=unparsed_args)
        except Exception as exc:
            raise self.SubcommandExecutionError(subcommand) from exc

    def _get_subcommand_fn(self, subcommand=None):
        return getattr(self, subcommand)

    def _call_subcommand_fn(self, subcommand_fn=None, parsed_args=None,
                            unparsed_args=None):
        return subcommand_fn(parsed_args=parsed_args,
                             unparsed_args=unparsed_args)
