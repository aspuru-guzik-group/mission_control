from mc.utils.mc_sandbox import McSandbox

from mc.flows.flow_engine import FlowEngine
from mc.runners.flow_runner import FlowRunner


def main():
    sandbox = McSandbox()
    mc_db = sandbox.mc_db
    flow_engine = FlowEngine()
    create_flows(mc_db=mc_db, flow_engine=flow_engine)
    flow_runner = FlowRunner(
        flow_engine=flow_engine,
        flow_record_client=sandbox.flow_record_client,
    )
    tick_stats = flow_runner.tick()
    while tick_stats['claimed'] > 0:
        tick_stats = flow_runner.tick()
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
        tasks.append({
            'key': task_id,
            'task_type': 'print',
            'task_params': {
                'msg': msg_tpl.format(task_id=task_id, flow_id=flow_id)
            }
        })
    flow_spec['tasks'] = tasks
    return flow_spec


if __name__ == '__main__':
    main()
