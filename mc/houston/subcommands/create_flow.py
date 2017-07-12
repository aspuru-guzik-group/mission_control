import json

from mc.flows.flow import Flow
from ._base_houston_subcommand import BaseHoustonSubcommand

class CreateFlowSubcommand(BaseHoustonSubcommand):
    def _run(self, args=None, kwargs=None, unparsed_args=None):
        flow_spec = {
            'tasks': [
                # *[
                    # {'task_type': 'print', 'task_params': {'msg': i}}
                    # for i in range(3)
                # ],
                # {
                    # 'task_type': 'tasks.example_countdown',
                    # 'task_params': {'countdown_start': 3}
                # },
                {
                    'task_type': 'job',
                    'task_params': {
                        'job_type': 'job_modules.example_echo',
                        'job_params': {
                            'message' :'Hello echo!'
                        }
                    }
                }
            ]
        }
        flow_dict = Flow.from_flow_spec(flow_spec=flow_spec).to_flow_dict()
        flow_record = self.utils.mc_dao.create_item(item_type='Flow',
                                                    item_kwargs=flow_dict)
        print("Created flow_record:\n", json.dumps(flow_record, indent=2))

Subcommand = CreateFlowSubcommand
