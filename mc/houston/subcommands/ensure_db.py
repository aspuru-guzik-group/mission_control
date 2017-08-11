from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def _run(self):
        self.utils.ensure_db()
        return {'msg': 'ensured dbs'}
