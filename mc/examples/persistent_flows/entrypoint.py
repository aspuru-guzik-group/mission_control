import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.daos.sa_dao import SaDao

_DIR = os.path.dirname(__file__)

def main():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_uri = 'sqlite:///{db_file}'.format(db_file=db_file)
    dao = SaDao(db_uri=db_uri)
    dao.create_tables()

    flow_spec = generate_flow_spec()
    flow_engine = FlowEngine()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)

    flow_record = dao.create_item(item_type='Flow', kwargs={
        'serialization': flow_engine.serialize_flow(flow=flow)
    })
    flow_key = flow_record['key']

    fetched_flow_record_1 = dao.get_item_by_key(item_type='Flow', key=flow_key)
    reconstituted_flow_1 = flow_engine.deserialize_flow(
        serialized_flow=fetched_flow_record_1['serialization'])
    flow_engine.run_flow(flow=reconstituted_flow_1, check=True)
    dao.patch_item(item_type='Flow', key=flow_key, patches={
        'serialization': flow_engine.serialize_flow(reconstituted_flow_1),
        'status': reconstituted_flow_1.status
    })

    fetched_flow_record_2 = dao.get_item_by_key(item_type='Flow', key=flow_key)
    reconstituted_flow_2 = flow_engine.deserialize_flow(
        serialized_flow=fetched_flow_record_2['serialization'])
    print(reconstituted_flow_2.status)

def generate_flow_spec():
    flow_spec = {
        'label': 'example_flow',
        'tasks': [
            {
                'key': 'task_1',
                'task_type': 'print',
                'task_params': {'message': 'I am task_1.'},
                'precursors': ['ROOT'],
            },
            {
                'key': 'task_2',
                'task_type': 'print',
                'task_params': {'message': 'I am task_2.'},
                'precursors': ['task_1'],
            },
            {
                'key': 'task_3',
                'task_type': 'print',
                'task_params': {'message': 'I am task_3.'},
                'precursors': ['task_2'],
            }
        ]
    }
    return flow_spec

if __name__ == '__main__': main()
