import json

from mc.flows.flow import Flow
from ._base_houston_subcommand import BaseHoustonSubcommand

class CreateFlowSubcommand(BaseHoustonSubcommand):
    def add_arguments(self, parser=None):
        parser.add_argument('--flow_spec', type=json.loads)

    def _run(self):
        flow_spec = self.kwargs['flow_spec']
        flow_dict = Flow.from_flow_spec(flow_spec=flow_spec).to_flow_dict()
        flow_record = self.utils.mc_dao.create_item(
            item_type='flow', item_kwargs=flow_dict)
        print(json.dumps(flow_record))

Subcommand = CreateFlowSubcommand
