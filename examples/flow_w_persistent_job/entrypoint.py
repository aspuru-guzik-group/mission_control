from mc.utils.mc_sandbox import McSandbox


def main():
    sandbox = McSandbox()
    flow_spec = generate_flow_spec()
    flow = sandbox.flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
    task_ctx = setup_task_ctx(mc_db=sandbox.mc_db)
    sandbox.flow_engine.run_flow(flow=flow, task_ctx=task_ctx)


def generate_flow_spec():
    flow_spec = {
        'label': 'example_flow',
        'tasks': [
            {
                'key': 'task_1',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_1.'},
            },
            {
                'key': 'task_2',
                'task_type': 'mc.tasks.job',
                'task_params': {
                    'job_type': 'some.job.type'
                },
            },
        ]
    }
    return flow_spec


def setup_task_ctx(mc_db=None):
    class MyJobRecordClient(object):
        def __init__(self, mc_db=None):
            self.mc_db = mc_db

        def create_job_record(self, *args, job_kwargs=None, **kwargs):
            job_kwargs = {**(job_kwargs or {}), 'data': {},
                          'status': 'PENDING'}
            job_record = self.mc_db.create_item(
                item_type='job', item_kwargs=job_kwargs)
            job_meta = {'key': job_record['key']}
            return job_meta

        def get_job_record(self, *args, job_meta=None, **kwargs):
            job_record = self.mc_db.get_item_by_key(
                item_type='job', key=job_meta['key'])
            # Normally the job would be run by an external runner,
            # but we're faking it here.
            tick_job_record(job_record=job_record)
            self.mc_db.patch_item(
                item_type='job', key=job_meta['key'],
                patches={'data': job_record['data'],
                         'status': job_record['status']}
            )
            return job_record

    def tick_job_record(job_record=None):
        job_record['data'].setdefault('tick_count', 0)
        job_record['data']['tick_count'] += 1
        if job_record['data']['tick_count'] < 2:
            print("Running job_record.")
            job_record['status'] = 'RUNNING'
        else:
            print("Completing job_record.")
            job_record['status'] = 'COMPLETED'

    task_ctx = {'mc.job_record_client': MyJobRecordClient(mc_db=mc_db)}
    return task_ctx


if __name__ == '__main__':
    main()
