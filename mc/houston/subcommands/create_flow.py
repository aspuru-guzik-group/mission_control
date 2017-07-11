def create_flow(self, args=None, kwargs=None, unparsed_args=None):
    mc_dao = self._get_mc_dao()
    flow_spec = {
        'tasks': [
            *[
                {'task_type': 'print', 'task_params': {'msg': i}}
                for i in range(3)
            ],
            {
                'task_type': 'tasks.example_countdown',
                'task_params': {'countdown_start': 3}
            },
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
    flow_record = mc_dao.create_item(item_type='Flow',
                                     item_kwargs=flow_dict)
    print("Created flow_record: ", flow_record)
