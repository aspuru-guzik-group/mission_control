from mc.utils.mc_sandbox import McSandbox


def main():
    sandbox = McSandbox()
    mc_dao = sandbox.mc_dao
    flow_engine = sandbox.flow_engine
    create_flows(mc_dao=mc_dao, flow_engine=flow_engine, n=3)

    queue_record = mc_dao.create_item(
        item_type='Queue',
        kwargs={'queue_spec': {'item_type': 'Flow'}
    })
    queue_key = queue_record['key']

    claimed_flow_records = claim_flow_records(mc_dao=mc_dao, queue_key=queue_key)
    while claimed_flow_records:
        for flow_record in claimed_flow_records:
            reconstituted_flow = flow_engine.flow_dict_to_flow(
                flow_dict=flow_record)
            flow_engine.tick_flow(flow=reconstituted_flow)
            mc_dao.patch_item(
                item_type='Flow',
                key=flow_record['key'],
                patches={
                    **flow_engine.flow_to_flow_dict(flow=reconstituted_flow),
                    'claimed': False
                }
            )
        claimed_flow_records = claim_flow_records(mc_dao=mc_dao,
                                                  queue_key=queue_key)
    print("No more flows to claimed.")

def create_flows(mc_dao=None, flow_engine=None, n=3):
    for i in range(n):
        flow_spec = generate_flow_spec(params={
            'flow_id': 'flow_%s' % i,
            'num_tasks': 3
        })
        flow = flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
        mc_dao.create_item(item_type='Flow',
                           kwargs=flow_engine.flow_to_flow_dict(flow=flow))

def generate_flow_spec(params=None):
    msg_tpl = "I am task '{task_id}' in flow '{flow_id}'"
    flow_id = params['flow_id']
    num_tasks = params['num_tasks']
    flow_spec = {'label': flow_id}
    tasks = []
    for i in range(num_tasks):
        task_id = 'task_%s' % i
        if i == 0: precursor = 'ROOT'
        else: precursor = tasks[i-1]['key']
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

def claim_flow_records(mc_dao=None, queue_key=None):
    return mc_dao.claim_queue_items(queue_key=queue_key)['items']

if __name__ == '__main__': main()
