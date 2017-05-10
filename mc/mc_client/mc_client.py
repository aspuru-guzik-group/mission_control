import logging


class MissionControlClient(object):
    def __init__(self, dao=None, logger=None):
        self.dao = dao
        self.logger = logger or logging

    def create_queue(self, queue=None):
        return self.dao.create_queue(queue=queue)

    def claim_queue_items(self, queue_key=None, params=None):
        return self.dao.claim_queue_items(queue_key=queue_key, params=params)

    def create_flow(self, flow=None):
        return self.dao.create_flow(flow=flow)

    def get_flows(self, query=None):
        return self.dao.get_flows(query=query)

    def get_flow_by_uuid(self, uuid=None):
        flows = self.get_flows(query={
            'filters': [
                {'field': 'uuid', 'operator': '=', 'value': uuid}
            ]
        })
        if len(flows) == 1: return flows[0]
        else: return None

    def patch_flows(self, keyed_patches=None):
        keyed_results = {}
        for key, patches_for_key in (keyed_patches or {}).items():
            try: result = self.patch_flow(key=key, patches=patches_for_key)
            except Exception as exc: result = {"exception", exc}
            keyed_results[key]= result
        return keyed_results

    def patch_flow(self, key=None, patches=None):
        return self.dao.patch_flow(key=key, patches=patches)

    def create_job(self, job=None):
        return self.dao.create_job(job=job)

    def get_jobs(self, query=None):
        return self.dao.get_jobs(query=query)

    def get_job_by_uuid(self, uuid=None):
        jobs = self.get_jobs(query={
            'filters': [
                {'field': 'uuid', 'operator': '=', 'value': uuid}
            ]
        })
        if len(jobs) == 1: return jobs[0]
        else: return None

    def patch_jobs(self, keyed_patches=None):
        keyed_results = {}
        for key, patches_for_key in (keyed_patches or {}).items():
            try: result = self.patch_job(key=key, patches=patches_for_key)
            except Exception as exc: result = {"exception", exc}
            keyed_results[key]= result
        return keyed_results

    def patch_job(self, key=None, patches=None):
        return self.dao.patch_job(key=key, patches=patches)

    def flush_mc_db(self):
        return self.dao.flush_mc_db()
