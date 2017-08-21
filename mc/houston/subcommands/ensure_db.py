from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def _run(self):
        self.houston.db.ensure_tables()
        return {'msg': 'ensured db'}
