import unittest

from .. import db


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.db = db.Db(db_uri='sqlite://', ensure=True)


class FlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_kwargs = {'label': 'initial_label'}

    def _create_flow(self):
        return self.db.create_item(
            item_type='flow', item_kwargs=self.flow_kwargs)

    def test_create_flow(self):
        created_flow = self._create_flow()
        self.assertEqual(created_flow['label'], self.flow_kwargs['label'])
        self.assertTrue(created_flow['key'] is not None)

    def test_patch_flow(self):
        created_flow = self._create_flow()
        patched_flow = self.db.patch_item(
            item_type='flow',
            key=created_flow['key'],
            patches={'label': 'label2'}
        )
        self.assertEqual(patched_flow['label'], 'label2')


class JobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_kwargs = {
            'label': 'initial_label',
            'data': {'some': 'data'},
            'job_type': 'some.job_type',
            'job_params': {'some': 'job_params'},
            'cfg': {'some': 'cfg'},
            'parent_key': 'some_parent_key',
        }

    def _create_job(self):
        return self.db.create_item(
            item_type='job', item_kwargs=self.job_kwargs)

    def test_create_job(self):
        created_job = self._create_job()
        for k, v in self.job_kwargs.items():
            self.assertEqual(created_job[k], v)
        self.assertTrue(created_job['key'] is not None)

    def test_patch_job(self):
        created_job = self._create_job()
        patches = {
            'label': 'label2',
            'data': {'some': 'data2'},
            'job_type': 'some.job_type2',
            'job_params': {'some': 'job_params2'},
            'cfg': {'some': 'cfg2'},
            'parent_key': 'some_parent_key2',
        }
        patched_job = self.db.patch_item(
            item_type='job', key=created_job['key'], patches=patches)
        for k, v in patches.items():
            self.assertEqual(patched_job[k], v)


class JobQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_kwargs = {'label': 'initial_label'}

    def _create_queue(self, queue_kwargs=None):
        queue_kwargs = queue_kwargs or self.queue_kwargs
        return self.db.create_item(
            item_type='queue', item_kwargs=queue_kwargs)

    def test_create_queue(self):
        created_queue = self._create_queue()
        self.assertEqual(created_queue['label'], self.queue_kwargs['label'])
        self.assertTrue(created_queue['key'] is not None)

    def test_claim_queue_items(self):
        queue_kwargs = {'label': 'initial_label',
                        'queue_spec': {'item_type': 'job'}}
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        flows = [
            self.db.create_item(
                item_type='job', item_kwargs={'label': 'flow_%s' % i})
            for i in range(3)
        ]
        claimed_flows = (
            self.db.claim_queue_items(queue_key=queue['key'])['items'])
        actual_keys = sorted([flow['key'] for flow in claimed_flows])
        expected_keys = sorted([flow['key'] for flow in flows])
        self.assertEqual(actual_keys, expected_keys)


class FlowQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_kwargs = {'label': 'initial_label'}

    def _create_queue(self, queue_kwargs=None):
        queue_kwargs = queue_kwargs or self.queue_kwargs
        return self.db.create_item(
            item_type='queue', item_kwargs=queue_kwargs)

    def test_create_queue(self):
        created_queue = self._create_queue()
        self.assertEqual(created_queue['label'], self.queue_kwargs['label'])
        self.assertTrue(created_queue['key'] is not None)

    def _claim_flows(self, queue=None):
        return self.db.claim_queue_items(queue_key=queue['key'])['items']

    def _release_flows(self, flows=None):
        self.db.patch_items(
            item_type='flow',
            keyed_patches={flow['key']: {'claimed': False} for flow in flows}
        )

    def test_claim_queue_items(self):
        queue_kwargs = {'label': 'initial_label',
                        'queue_spec': {'item_type': 'flow'}}
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        flows = [
            self.db.create_item(
                item_type='flow', item_kwargs={'label': 'flow_%s' % i})
            for i in range(3)
        ]
        claimed_flows = self._claim_flows(queue=queue)
        actual_keys = sorted([flow['key'] for flow in claimed_flows])
        expected_keys = sorted([flow['key'] for flow in flows])
        self.assertEqual(actual_keys, expected_keys)

    def test_locks(self):
        queue_kwargs = {'label': 'initial_label',
                        'queue_spec': {'item_type': 'flow'}}
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        unlocked_flow = self.db.create_item(
            item_type='flow', item_kwargs={'label': 'unlocked'})
        flow_to_unlock = self.db.create_item(
            item_type='flow',
            item_kwargs={'label': 'to_unlock',
                         'num_tickable_tasks': 1}
        )
        locker_key = 'the_locker'
        self.db.create_lock(
            lockee_key=flow_to_unlock['key'], locker_key=locker_key)
        claimed_flows = self._claim_flows(queue=queue)
        self.assert_flow_lists_match(claimed_flows, [unlocked_flow])
        self._release_flows(flows=claimed_flows)
        self.db.release_locks(locker_keys=[locker_key])
        claimed_flows = self._claim_flows(queue=queue)
        self.assert_flow_lists_match(
            claimed_flows, [unlocked_flow, flow_to_unlock])
        self._release_flows(flows=claimed_flows)

    def assert_flow_lists_match(self, flows_a, flows_b):
        self.assertEqual(sorted([flow['key'] for flow in flows_a]),
                         sorted([flow['key'] for flow in flows_b]))
