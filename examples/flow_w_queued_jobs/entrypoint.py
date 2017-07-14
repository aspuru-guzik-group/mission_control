from mc.utils.mc_sandbox import McSandbox


def main():
    sandbox = McSandbox()
    flow_spec = generate_flow_spec()
    sandbox.flow_record_client.create_flow_record_from_flow_spec(
        flow_spec=flow_spec)
    while sandbox.has_incomplete_items():
        sandbox.flow_runner.tick()
        claimed_jobs = sandbox.mc_dao.claim_queue_items(
            queue_key=sandbox.queues['job']['key'])['items']
        for job in claimed_jobs: run_job(job=job, mc_dao=sandbox.mc_dao)
    print("No more flows to claim.")

def generate_flow_spec():
    flow_spec = {
        'label': 'example_flow',
        'tasks': []
    }
    for i in range(3):
        job_type = 'job.%s' % i
        flow_spec['tasks'].append({
            'key': 'job_{job_type}_task'.format(job_type=job_type),
            'task_type': 'mc.tasks.job',
            'task_params': {'job_type': job_type},
            'precursors': ['ROOT'],
        })
    return flow_spec

def run_job(job=None, mc_dao=None):
    print("running job '{job_type}'".format(job_type=job['job_type']))
    job['status'] = 'COMPLETED'
    mc_dao.patch_item(item_type='job', key=job['key'], patches={
        'data': job['data'],
        'status': job['status'],
    })

if __name__ == '__main__': main()
