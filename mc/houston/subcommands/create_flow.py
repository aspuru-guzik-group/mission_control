import json
import os

from mc.flows.flow import Flow
from ._base_houston_subcommand import BaseHoustonSubcommand

class CreateFlowSubcommand(BaseHoustonSubcommand):
    def add_arguments(self, parser=None):
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--flow_spec', type=json.loads)
        group.add_argument('--flow_spec_file', type=self._json_path)

    def _json_path(self, path):
        with open(os.path.expanduser(path)) as f: return json.load(f)

    def _run(self):
        flow_spec = self.parsed_args.get('flow_spec')
        if flow_spec is None: self.parsed_args.get('flow_spec_file')
        flow_dict = Flow.from_flow_spec(flow_spec=flow_spec).to_flow_dict()
        flow_record = self.utils.mc_dao.create_item(
            item_type='flow', item_kwargs=flow_dict)
        print(json.dumps(flow_record))

Subcommand = CreateFlowSubcommand
