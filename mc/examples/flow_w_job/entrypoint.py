import os

from mc.flow_engines.flow_engine import FlowEngine

_DIR = os.path.dirname(__file__)

def main():
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
    flow_engine = FlowEngine()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)
    task_ctx = setup_task_ctx()
    flow_engine.run_flow(flow=flow, task_ctx=task_ctx)

def setup_task_ctx():
    def create_job(*args, **kwargs):
        print("creating job, args: {args}, kwargs: {kwargs}".format(
            args=args, kwargs=kwargs))
        job_meta = {'key': 'some_key'}
        return job_meta

    def get_job(*args, job_meta=None, **kwargs):
        print("getting job, job_meta: {job_meta}".format(job_meta=job_meta))
        fake_job = {
            'status': 'COMPLETED',
            'data': {
                'artifact': 'fake artifact',
                'std_logs': 'fake_std_logs'
            }
        }
        return fake_job

    task_ctx = {
        'mc.tasks.job.create_job': create_job,
        'mc.tasks.job.get_job': get_job,
    }
    return task_ctx

if __name__ == '__main__': main()
