class JobRecordClient(object):
    def __init__(self, mc_dao=None, queue_key=None, use_locks=False):
        self.mc_dao = mc_dao
        self.queue_key = queue_key
        self.use_locks = use_locks

    def create_job_record(self, *args, job_kwargs=None, **kwargs):
        job_kwargs = {**(job_kwargs or {}), 'status': 'PENDING'}
        job_record = self.mc_dao.create_item(item_type='Job', kwargs=job_kwargs)
        if self.use_locks:
            parent_key = job_kwargs.get('data', {}).get('parent_key')
            if parent_key:
                self.mc_dao.create_lock(lockee_key=parent_key,
                                        locker_key=job_record['key'])
        job_meta = {'key': job_record['key']}
        return job_meta

    def get_job_record(self, job_meta=None):
        return self.mc_dao.get_item_by_key(item_type='Job', key=job_meta['key'])

    def claim_job_records(self, *args, **kwargs):
        return self.mc_dao.claim_queue_items(queue_key=self.queue_key)['items']
    
    def patch_job_records(self, keyed_patches=None):
        patched = self.mc_dao.patch_items(item_type='Job',
                                          keyed_patches=keyed_patches)
        if self.use_locks:
            locker_keys = [job_key for job_key, patches in keyed_patches.items()
                           if patches.get('status') in {'FAILED', 'COMPLETED'}]
            if locker_keys: self.mc_dao.release_locks(locker_keys=locker_keys)
        return patched
