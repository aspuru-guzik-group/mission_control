from ._base_houston_subcommand import BaseHoustonSubcommand


class SanityCheckSubcommand(BaseHoustonSubcommand):
    def _run(self): print("sane!")

Subcommand = SanityCheckSubcommand
