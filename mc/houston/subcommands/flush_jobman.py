from ._base_houston_subcommand import BaseHoustonSubcommand

class FlushJobmanSubcommand(BaseHoustonSubcommand):
    def _run(self): self.utils.jobman.flush()

Subcommand = FlushJobmanSubcommand
