import json

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Query MC items"

    def add_arguments(self, parser=None):
        defaults = self._get_defaults()
        parser.add_argument('--item_type')
        parser.add_argument(
            '--query',
            help=('JSON string representing query'),
            type=json.loads,
            default=defaults['query']
        )

    def _get_defaults(self):
        return {
            'query': {}
        }

    def _run(self):
        return self.houston.utils.db.delete_items(
            item_type=self.parsed_args['item_type'],
            query=self.parsed_args['query']
        )
