import json

from ._base_houston_subcommand import BaseHoustonSubcommand


class InfoSubcommand(BaseHoustonSubcommand):
    def add_arguments(self, parser=None):
        pass

    def _run(self, args=None, kwargs=None, unparsed_args=None):
        info = self._get_infos_for_mc_record_types()
        print(json.dumps(info))

    def _get_infos_for_mc_record_types(self):
        infos_for_types = {
            record_type: self._get_info_for_mc_record_type(
                record_type=record_type)
            for record_type in ['Flow', 'Job']
        }
        return infos_for_types

    def _get_info_for_mc_record_type(self, record_type=None):
        info_for_type = {
            'count': len(self.utils.mc_dao.get_items(item_type=record_type))
        }
        return info_for_type


Subcommand = InfoSubcommand
