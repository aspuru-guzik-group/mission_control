import json

from ._base_houston_subcommand import BaseHoustonSubcommand


class InfoSubcommand(BaseHoustonSubcommand):
    def add_arguments(self, parser=None):
        parser.add_argument('--key', help="key for a single object")
        parser.add_argument('--indent', type=int, help="indent output")

    def _run(self):
        key = self.parsed_args.get('key')
        if key: info = self._get_info_for_key(key=key)
        else: info = self._get_mc_record_type_summaries()
        print(json.dumps(info, indent=self.parsed_args.get('indent')))

    def _get_info_for_key(self, key=None):
        record_type = key.split(':')[0]
        if record_type not in ['flow', 'job', 'queue', 'lock']:
            raise Exception("Invalid key '%s'" % key)
        return self.utils.mc_dao.get_item_by_key(item_type=record_type,
                                                 key=key)

    def _get_mc_record_type_summaries(self):
        summaries = {
            record_type: self._get_mc_record_type_summary(
                record_type=record_type)
            for record_type in ['flow', 'job']
        }
        return summaries

    def _get_mc_record_type_summary(self, record_type=None):
        summary = {
            'count': len(self.utils.mc_dao.get_items(item_type=record_type))
        }
        return summary


Subcommand = InfoSubcommand
