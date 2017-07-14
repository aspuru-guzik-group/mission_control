import json

from ._base_houston_subcommand import BaseHoustonSubcommand


class EnsureDbSubcommand(BaseHoustonSubcommand):
    def _run(self):
        self.utils.ensure_db()
        print(json.dumps({'msg': 'ensured dbs'}))


Subcommand = EnsureDbSubcommand
