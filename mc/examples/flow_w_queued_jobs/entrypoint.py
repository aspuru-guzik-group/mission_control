import json
import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.mc_daos.sa_dao import SaDao
from mc.runners.flow_runner import FlowRunner

_DIR = os.path.dirname(__file__)

def main():
    dao = setup_dao()
    flow_engine = FlowEngine()
    create_flow(dao=dao, flow_engine=flow_engine)

    flow_queue_key = dao.create_item(item_type='Queue', kwargs={
        'queue_spec': json.dumps({'item_type': 'Flow'})
    })['key']
    job_queue_key = dao.create_item(item_type='Queue', kwargs={
        'queue_spec': json.dumps({'item_type': 'Job'})
    })['key']

    flow_runner = FlowRunner(
        flow_client=FlowRunner.generate_flow_client_from_mc_dao(
            dao=dao, queue_key=flow_queue_key),
        task_context=setup_task_context(dao=dao)
    )

    tick_stats = flow_runner.tick()
    while tick_stats['claimed'] > 0:
        tick_stats = flow_runner.tick()
        claimed_job_records = dao.claim_queue_items(
            queue_key=job_queue_key)['items']
        for job_record in claimed_job_records:
            run_job_record(job_record=job_record, dao=dao)
    print("No more flows to claim.")

def setup_dao():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_uri = 'sqlite:///{db_file}'.format(db_file=db_file)
    dao = SaDao(db_uri=db_uri)
    dao.create_tables()
    return dao

def create_flow(dao=None, flow_engine=None):
    flow_spec = generate_flow_spec()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)
    dao.create_item(item_type='Flow', kwargs={
        'serialization': flow_engine.serialize_flow(flow=flow)
    })

def generate_flow_spec():
    flow_spec = {
        'label': 'example_flow',
        'task_specs': []
    }
    for i in range(3):
        job_id = 'job_%s' % i
        flow_spec['task_specs'].append({
            'task' : {
                'key': 'job_{job_id}_task'.format(job_id=job_id),
                'task_type': 'mc.tasks.job',
                'task_params': {'job_spec': {'job_id': job_id}},
            },
            'precursors': ['ROOT'],
        })
    return flow_spec

def setup_task_context(dao=None):
    def create_job(*args, job_kwargs=None, **kwargs):
        job = {**job_kwargs, 'data': {}, 'status': 'PENDING'}
        job_record = dao.create_item(item_type='Job', kwargs={
            'serialization': dao.serialize_value(job),
            'status': job['status'],
        })
        job_meta = {'key': job_record['key']}
        return job_meta

    def get_job(*args, job_meta=None, **kwargs):
        job_record = dao.get_item_by_key(item_type='Job', key=job_meta['key'])
        job = dao.deserialize_value(job_record['serialization'])
        return job

    task_context = {
        'mc.tasks.job.create_job': create_job,
        'mc.tasks.job.get_job': get_job,
    }
    return task_context

def run_job_record(job_record=None, dao=None):
    job = dao.deserialize_value(job_record['serialization'])
    print("running job '{job_id}'".format(job_id=job['job_spec']['job_id']))
    job['status'] = 'COMPLETED'
    dao.patch_item(item_type='Job', key=job_record['key'], patches={
        'serialization': dao.serialize_value(job),
        'status': job['status']
    })

if __name__ == '__main__': main()
