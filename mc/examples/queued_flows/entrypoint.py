import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.daos.sa_dao import SaDao

_DIR = os.path.dirname(__file__)

def main():
    dao = setup_dao()
    flow_engine = FlowEngine()
    create_flows(dao=dao, flow_engine=flow_engine, n=3)

    queue_record = dao.create_item(item_type='Queue', kwargs={
        'queue_spec': {'item_type': 'Flow'}
    })
    queue_key = queue_record['key']

    claimed_flow_records = claim_flow_records(dao=dao, queue_key=queue_key)
    while claimed_flow_records:
        for flow_record in claimed_flow_records:
            reconstituted_flow = flow_engine.deserialize_flow(
                serialized_flow=flow_record['serialization'])
            flow_engine.tick_flow(flow=reconstituted_flow)
            dao.patch_item(
                item_type='Flow',
                key=flow_record['key'],
                patches={
                    'serialization': flow_engine.serialize_flow(
                        reconstituted_flow),
                    'claimed': False,
                    'status': reconstituted_flow.status
                }
            )
        claimed_flow_records = claim_flow_records(dao=dao, queue_key=queue_key)
    print("No more flows to claimed.")

def setup_dao():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_uri = 'sqlite:///{db_file}'.format(db_file=db_file)
    dao = SaDao(db_uri=db_uri)
    dao.create_tables()
    return dao

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

def claim_flow_records(dao=None, queue_key=None):
    return dao.claim_queue_items(queue_key=queue_key)['items']

if __name__ == '__main__': main()
