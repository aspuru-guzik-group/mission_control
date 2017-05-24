import os
import tempfile

from mc.flow_engines.flow_engine import FlowEngine
from mc.mc_daos.sa_dao import SaDao

_DIR = os.path.dirname(__file__)

def main():
    dao = setup_dao()
    flow_engine = FlowEngine()
    flow_spec = generate_flow_spec()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)
    task_context = setup_task_context(dao=dao)
    flow_engine.run_flow(flow=flow, task_context=task_context)

def setup_dao():
    db_file = tempfile.mkstemp(suffix='.sqlite.db')[1]
    db_uri = 'sqlite:///{db_file}'.format(db_file=db_file)
    dao = SaDao(db_uri=db_uri)
    dao.create_tables()
    return dao

def generate_flow_spec():
    flow_spec = {
        'label': 'example_flow',
        'task_specs': [
            {
                'task' : {
                    'key': 'task_1',
                    'task_type': 'print',
                    'task_params': {'message': 'I am task_1.'},
                },
                'precursors': ['ROOT'],
            },
            {
                'task' : {
                    'key': 'task_2',
                    'task_type': 'mc.tasks.job',
                    'task_params': {
                        'job_spec': 'some job spec'
                    },
                },
                'precursors': ['task_1'],
            },
        ]
    }
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
        # Normally the job would be run by an external runner, but we're faking
        # it here.
        tick_job(job=job)
        dao.patch_item(item_type='Job', key=job_meta['key'], patches={
            'serialization': dao.serialize_value(job),
            'status': job['status']
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

    task_context = {
        'mc.tasks.job.create_job': create_job,
        'mc.tasks.job.get_job': get_job,
    }
    return task_context

if __name__ == '__main__': main()
