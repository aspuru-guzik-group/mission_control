import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.mc_daos.sa_dao import SaDao
from mc.runners.flow_runner import FlowRunner

_DIR = os.path.dirname(__file__)

def main():
    dao = setup_dao()
    flow_engine = FlowEngine()
    queue_record = create_queue(dao=dao)
    create_flows(dao=dao, flow_engine=flow_engine)
    flow_runner = FlowRunner(
        flow_engine=flow_engine,
        flow_client=FlowRunner.generate_flow_client_from_mc_dao(
            dao=dao, queue_key=queue_record['key'])
    )
    tick_stats = flow_runner.tick()
    while tick_stats['claimed'] > 0: tick_stats = flow_runner.tick()
    print("No more flows to claimed.")

def setup_dao():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_uri = 'sqlite:///{db_file}'.format(db_file=db_file)
    dao = SaDao(db_uri=db_uri)
    dao.create_tables()
    return dao

def create_queue(dao=None):
    queue_spec = {'item_type': 'Flow'}
    queue_record = dao.create_item(item_type='Queue', kwargs={
        'queue_spec': dao.serialize_value(queue_spec)
    })
    return queue_record

def create_flows(dao=None, flow_engine=None, n=3):
    for i in range(n):
        flow_spec = generate_flow_spec(params={
            'flow_id': 'flow_%s' % i,
            'num_tasks': 3
        })
        flow = flow_engine.generate_flow(flow_spec=flow_spec)
        dao.create_item(item_type='Flow', kwargs={
            'serialization': flow_engine.serialize_flow(flow=flow)
        })

def generate_flow_spec(params=None):
    msg_tpl = "I am task '{task_id}' in flow '{flow_id}'"
    flow_id = params['flow_id']
    num_tasks = params['num_tasks']
    flow_spec = {'label': flow_id}
    task_specs = []
    for i in range(num_tasks):
        task_id = 'task_%s' % i
        if i == 0: precursor = 'ROOT'
        else: precursor = task_specs[i-1]['task']['key']
        task_specs.append({
            'task' : {
                'key': task_id,
                'task_type': 'print',
                'task_params': {
                    'message': msg_tpl.format(task_id=task_id, flow_id=flow_id)
                },
            },
            'precursors': [precursor]
        })
    flow_spec['task_specs'] = task_specs
    return flow_spec

if __name__ == '__main__': main()
