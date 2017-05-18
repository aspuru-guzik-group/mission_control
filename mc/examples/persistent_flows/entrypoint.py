import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.mc_client.mc_client import MissionControlClient
from mc.mc_client.dao.dj_dao import DjDao

_DIR = os.path.dirname(__file__)

def main():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_params = {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': db_file
    }
    mc_client = MissionControlClient(dao=DjDao())
    flow_spec = generate_flow_spec()
    flow_engine = FlowEngine()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)
    flow_record = mc_client.create_flow(flow_kwargs={
        'serialization': flow_engine.serialize_flow(flow=flow)
    })
    flow_key = flow_record['uuid']

    fetched_flow_record_1 = mc_client.get_flow(flow_key=flow_key)
    reconstituted_flow_1 = flow_engine.deserialize_flow(
        serialized_flow=fetched_flow_record_1['serialization'])
    flow_engine.run_flow(flow=reconstituted_flow_1, check=True)
    mc_client.patch_flow(
        flow_key=flow_key,
        patches={
            'serialization': flow_engine.serialize_flow(reconstituted_flow_1),
            'status': reconstituted_flow_1['status']
        }
    )

    fetched_flow_record_2 = mc_client.get_flow(flow_key=flow_key)
    reconstituted_flow_2 = flow_engine.deserialize_flow(
        serialized_flow=fetched_flow_record_2['serialization'])
    fetched_flow_record_2 = mc_client.get_flow(flow_key=flow_key)
    print(reconstituted_flow_2.status)

def generate_mc_client(mc_db_params=None):
    pass

def generate_flow_spec():
    pass

if __name__ == '__main__': main()
