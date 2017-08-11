import json

from mc.flows.flow import Flow
from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Create a flow"

    def add_arguments(self, parser):
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--flow_spec', type=json.loads)
        group.add_argument('--flow_spec_file', type=self.json_path_arg)

    def _get_defaults(self):
        return {
            'flow_spec': {}
        }

    def _ensure_parsed_args(self, *args, **kwargs):
        super()._ensure_parsed_args(*args, **kwargs)
        flow_spec = self.parsed_args.get('flow_spec')
        if flow_spec is None:
            flow_spec = self.parsed_args.get('flow_spec_file')

    def _run(self):
        flow_spec = self.parsed_args['flow_spec']
        flow_dict = Flow.from_flow_spec(flow_spec=flow_spec).to_flow_dict()
        flow_record = self.utils.db.create_item(
            item_type='flow', item_kwargs=flow_dict)
        return flow_record
