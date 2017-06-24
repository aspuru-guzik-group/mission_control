from mc.flows.flow_engine import FlowEngine

class FlowRecordClient(object):
    def __init__(self, mc_dao=None, queue_key=None, use_locks=False):
        self.mc_dao = mc_dao
        self.queue_key = queue_key
        self.use_locks = use_locks

    def create_flow_record(self, flow_kwargs=None):
        flow_record = self.mc_dao.create_item(item_type='Flow',
                                              kwargs=flow_kwargs)
        if self.use_locks:
            parent_key = flow_kwargs.get('data', {}).get('parent_key')
            if parent_key: self.mc_dao.create_lock(
                lockee_key=parent_key, locker_key=flow_record['key'])
        flow_meta = {'key': flow_record['key']}
        return flow_meta

    def get_flow_record(self, flow_meta=None):
        return self.mc_dao.get_item_by_key(item_type='Flow',
                                           key=flow_meta['key'])

    def create_flow_record_from_flow_spec(self, flow_spec=None):
        flow = FlowEngine.flow_spec_to_flow(flow_spec=flow_spec)
        flow_dict = FlowEngine.flow_to_flow_dict(flow=flow)
        return self.create_flow_record(flow_kwargs=flow_dict)

    def claim_flow_records(self):
        claimed = self.mc_dao.claim_queue_items(
            queue_key=self.queue_key)['items']
        return claimed

    def patch_and_release_flow_record(self, flow_record=None, patches=None):
        patched =  self.mc_dao.patch_item(item_type='Flow',
                                          key=flow_record['key'],
                                          patches={'claimed': False, **patches})
        if self.use_locks:
            if patches.get('status') in {'FAILED', 'COMPLETED'}:
                self.mc_dao.release_locks(locker_keys=[flow_record['key']])
        return patched
