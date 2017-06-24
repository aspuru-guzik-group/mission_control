from mc.utils.mc_sandbox import McSandbox


def main():
    sandbox = McSandbox()
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
    flow = sandbox.flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
    task_ctx = setup_task_ctx()
    sandbox.flow_engine.run_flow(flow=flow, task_ctx=task_ctx)

def setup_task_ctx():
    class MyJobRecordClient(object):
        def create_job_record(self, *args, **kwargs):
            print("creating job_record")
            job_meta = {'key': 'some_key'}
            return job_meta

        def get_job_record(self, *args, job_meta=None, **kwargs):
            print("getting job_record,"
                  " job_meta: {job_meta}".format(job_meta=job_meta))
            fake_job = {
                'status': 'COMPLETED',
                'data': {
                    'artifact': 'fake artifact',
                    'std_logs': 'fake_std_logs'
                }
            }
            return fake_job

    task_ctx = {'mc.job_record_client': MyJobRecordClient()}
    return task_ctx

if __name__ == '__main__': main()
