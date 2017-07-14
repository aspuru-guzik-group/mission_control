from mc.utils.mc_sandbox import McSandbox

def main():
    sandbox = McSandbox()
    flow_engine = sandbox.flow_engine

    flow_spec = generate_flow_spec()
    flow = flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
    flow_record = sandbox.mc_dao.create_item(
        item_type='flow',
        item_kwargs=flow_engine.flow_to_flow_dict(flow=flow)
    )
    flow_key = flow_record['key']

    fetched_flow_record_1 = sandbox.mc_dao.get_item_by_key(item_type='flow',
                                                           key=flow_key)
    reconstituted_flow_1 = flow_engine.flow_dict_to_flow(
        flow_dict=fetched_flow_record_1)
    flow_engine.run_flow(flow=reconstituted_flow_1)
    sandbox.mc_dao.patch_item(
        item_type='flow', key=flow_key,
        patches=flow_engine.flow_to_flow_dict(flow=reconstituted_flow_1)
    )

    fetched_flow_record_2 = sandbox.mc_dao.get_item_by_key(item_type='flow',
                                                           key=flow_key)
    reconstituted_flow_2 = flow_engine.flow_dict_to_flow(
        flow_dict=fetched_flow_record_2)
    print(reconstituted_flow_2.status)

def generate_flow_spec():
    flow_spec = {
        'label': 'example_flow',
        'tasks': [
            {
                'key': 'task_1',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_1.'},
                'precursors': ['ROOT'],
            },
            {
                'key': 'task_2',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_2.'},
                'precursors': ['task_1'],
            },
            {
                'key': 'task_3',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_3.'},
                'precursors': ['task_2'],
            }
        ]
    }
    return flow_spec

if __name__ == '__main__': main()
