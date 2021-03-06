from mc.utils.mc_sandbox import McSandbox


def main():
    sandbox = McSandbox()
    mc_db = sandbox.mc_db
    flow_engine = sandbox.flow_engine
    create_flows(mc_db=mc_db, flow_engine=flow_engine, n=3)

    queue_record = mc_db.create_item(
        item_type='queue',
        item_kwargs={'queue_spec': {'item_type': 'flow'}}
    )
    queue_key = queue_record['key']

    claimed_flow_records = claim_flow_records(mc_db=mc_db, queue_key=queue_key)
    while claimed_flow_records:
        for flow_record in claimed_flow_records:
            reconstituted_flow = flow_engine.flow_dict_to_flow(
                flow_dict=flow_record)
            flow_engine.tick_flow(flow=reconstituted_flow)
            mc_db.patch_item(
                item_type='flow',
                key=flow_record['key'],
                patches={
                    **flow_engine.flow_to_flow_dict(flow=reconstituted_flow),
                    'claimed': False
                }
            )
        claimed_flow_records = claim_flow_records(mc_db=mc_db,
                                                  queue_key=queue_key)
    print("No more flows to claimed.")


def create_flows(mc_db=None, flow_engine=None, n=3):
    for i in range(n):
        flow_spec = generate_flow_spec(params={
            'flow_id': 'flow_%s' % i,
            'num_tasks': 3
        })
        flow = flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
        mc_db.create_item(
            item_type='flow',
            item_kwargs=flow_engine.flow_to_flow_dict(flow=flow)
        )


def generate_flow_spec(params=None):
    msg_tpl = "I am task '{task_id}' in flow '{flow_id}'"
    flow_id = params['flow_id']
    num_tasks = params['num_tasks']
    flow_spec = {'label': flow_id}
    tasks = []
    for i in range(num_tasks):
        task_id = 'task_%s' % i
        precursor = 'ROOT' if (i == 0) else tasks[i-1]['key']
        tasks.append({
            'key': task_id,
            'task_type': 'print',
            'task_params': {
                'msg': msg_tpl.format(task_id=task_id, flow_id=flow_id)
            },
            'precursors': [precursor]
        })
    flow_spec['tasks'] = tasks
    return flow_spec


def claim_flow_records(mc_db=None, queue_key=None):
    return mc_db.claim_queue_items(queue_key=queue_key)['items']


if __name__ == '__main__':
    main()
