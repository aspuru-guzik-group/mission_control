import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.daos.sa_dao import SaDao

_DIR = os.path.dirname(__file__)

def main():
    dao = setup_dao()
    flow_engine = FlowEngine()
    flow_spec = generate_flow_spec()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)
    task_ctx = setup_task_ctx(dao=dao)
    flow_engine.run_flow(flow=flow, task_ctx=task_ctx)

def setup_dao():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_uri = 'sqlite:///{db_file}'.format(db_file=db_file)
    dao = SaDao(db_uri=db_uri)
    dao.create_tables()
    return dao

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
                'task_type': 'mc.tasks.job',
                'task_params': {
                    'job_spec': 'some job spec'
                },
                'precursors': ['task_1'],
            },
        ]
    }
    return flow_spec

def setup_task_ctx(dao=None):
    def create_job(*args, job_kwargs=None, **kwargs):
        job_kwargs = {**(job_kwargs or {}), 'data': {}, 'status': 'PENDING'}
        job_record = dao.create_item(item_type='Job', kwargs=job_kwargs)
        job_meta = {'key': job_record['key']}
        return job_meta

    def get_job(*args, job_meta=None, **kwargs):
        job = dao.get_item_by_key(item_type='Job', key=job_meta['key'])
        # Normally the job would be run by an external runner, but we're faking
        # it here.
        tick_job(job=job)
        dao.patch_item(item_type='Job', key=job_meta['key'], patches={
            'data': job['data'],
            'status': job['status'],
        })
        return job

    def tick_job(job=None):
        job['data'].setdefault('tick_count', 0)
        job['data']['tick_count'] += 1
        if job['data']['tick_count'] < 2:
            print("Running job.")
            job['status'] = 'RUNNING'
        else:
            print("Completing job.")
            job['status'] = 'COMPLETED'

    task_ctx = {
        'mc.tasks.job.create_job': create_job,
        'mc.tasks.job.get_job': get_job,
    }
    return task_ctx

if __name__ == '__main__': main()
