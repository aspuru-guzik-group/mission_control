from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def _run(self): print("sane!")
