from mc.flows.flow_engine import FlowEngine


class FlowRecordClient(object):
    def __init__(self, mc_db=None, queue_key=None, use_locks=False):
        """
        Args:
            mc_db (mc_db): mc_db.
            queue_key (str): key for the flow queue to use.
            use_locks (bool): if True, create a lock record when creating.
                a job_record. This can be used to add a lock reference to
                parent flows which create nested flows.
        """
        self.mc_db = mc_db
        self.queue_key = queue_key
        self.use_locks = use_locks

    def create_flow_record(self, flow_kwargs=None):
        """
        Args:
            flow_kwargs (dict): a dictionary of flow parameters.

        Returns:
            flow_meta (dict): a dictionary of metadata that can be passed to
                get_flow_record to retrieve a flow.
        """
        flow_record = self.mc_db.create_item(
            item_type='flow', item_kwargs=flow_kwargs)
        if self.use_locks:
            parent_key = flow_kwargs.get('parent_key')
            if parent_key:
                self.mc_db.create_lock(
                    lockee_key=parent_key, locker_key=flow_record['key'])
        flow_meta = {'key': flow_record['key']}
        return flow_meta

    def get_flow_record(self, flow_meta=None):
        """
        Args:
            flow_meta (dict): dictionary of metadata, as per
                create_flow_record.

        Returns:
            flow_record (dict): a flow_record.
        """
        return self.mc_db.get_item_by_key(
            item_type='flow', key=flow_meta['key'])

    def create_flow_record_from_flow_spec(self, flow_spec=None):
        """Create a flow record from flow spec.

        Converts flow_spec to flow via flow_engine.

        Args:
            flow_spec (dict): a flow_spec

        Returns:
            flow_record (dict): a flow_record.
        """
        flow = FlowEngine.flow_spec_to_flow(flow_spec=flow_spec)
        flow_dict = FlowEngine.flow_to_flow_dict(flow=flow)
        return self.create_flow_record(flow_kwargs=flow_dict)

    def claim_flow_records(self):
        """
        Returns:
            flow_records (dict): a list of flow_records.
        """
        claimed = self.mc_db.claim_queue_items(
            queue_key=self.queue_key)['items']
        return claimed

    def patch_and_release_flow_record(self, flow_record=None, patches=None):
        """Patches a flow_record and clears its 'claimed' value.

        Args:
            flow_record (dict): a partial flow_record dict containing at least
                the 'key' field.
            patches (dict): a dict of flow_record kwargs to update.

        Returns:
            flow_record (dict): a patched flow dict.
        """
        patched = self.mc_db.patch_item(
            item_type='flow', key=flow_record['key'],
            patches={'claimed': False, **patches}
        )
        if self.use_locks:
            if patches.get('status') in {'FAILED', 'COMPLETED'}:
                self.mc_db.release_locks(locker_keys=[flow_record['key']])
        return patched
