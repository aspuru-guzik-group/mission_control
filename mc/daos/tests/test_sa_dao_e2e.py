import unittest

from .. import sa_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        #from sqlalchemy import create_engine
        #engine = create_engine('sqlite:///:memory:', echo=True)
        self.dao = sa_dao.SaDao(db_uri='sqlite:///:memory:')
        self.dao.create_tables()

class FlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_kwargs = {'label': 'initial_label'}

    def _create_flow(self):
        return self.dao.create_item(item_type='Flow', kwargs=self.flow_kwargs)

    def test_create_flow(self):
        created_flow = self._create_flow()
        self.assertEqual(created_flow['label'], self.flow_kwargs['label'])
        self.assertTrue(created_flow['key'] is not None)

    def test_patch_flow(self):
        created_flow = self._create_flow()
        patched_flow = self.dao.patch_item(
            item_type='Flow',
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
            'job_spec': {'some': 'job_spec'}
        }

    def _create_job(self):
        return self.dao.create_item(item_type='Job', kwargs=self.job_kwargs)

    def test_create_job(self):
        created_job = self._create_job()
        for k, v in self.job_kwargs.items(): self.assertEqual(created_job[k], v)
        self.assertTrue(created_job['key'] is not None)

    def test_patch_job(self):
        created_job = self._create_job()
        patches = {
            'label': 'label2',
            'data': {'some': 'data2'},
            'job_spec': {'some': 'job_spec2'},
        }
        patched_job = self.dao.patch_item(
            item_type='Job', key=created_job['key'], patches=patches)
        for k, v in patches.items(): self.assertEqual(patched_job[k], v)

class JobQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_kwargs = {'label': 'initial_label'}

    def _create_queue(self, queue_kwargs=None):
        queue_kwargs = queue_kwargs or self.queue_kwargs
        return self.dao.create_item(item_type='Queue', kwargs=queue_kwargs)

    def test_create_queue(self):
        created_queue = self._create_queue()
        self.assertEqual(created_queue['label'], self.queue_kwargs['label'])
        self.assertTrue(created_queue['key'] is not None)

    def test_claim_queue_items(self):
        queue_kwargs = {'label': 'initial_label',
                        'queue_spec': {'item_type': 'Job'}}
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        flows = [
            self.dao.create_item(
                item_type='Job',
                kwargs={'label': 'flow_%s' % i}
            )
            for i in range(3)
        ]
        claimed_flows = \
                self.dao.claim_queue_items(queue_key=queue['key'])['items']
        actual_keys = sorted([flow['key'] for flow in claimed_flows])
        expected_keys = sorted([flow['key'] for flow in flows])
        self.assertEqual(actual_keys, expected_keys)

class FlowQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_kwargs = {'label': 'initial_label'}

    def _create_queue(self, queue_kwargs=None):
        queue_kwargs = queue_kwargs or self.queue_kwargs
        return self.dao.create_item(item_type='Queue', kwargs=queue_kwargs)

    def test_create_queue(self):
        created_queue = self._create_queue()
        self.assertEqual(created_queue['label'], self.queue_kwargs['label'])
        self.assertTrue(created_queue['key'] is not None)

    def _claim_flows(self, queue=None):
        return self.dao.claim_queue_items(queue_key=queue['key'])['items']

    def _release_flows(self, flows=None):
        self.dao.patch_items(
            item_type='Flow',
            keyed_patches={flow['key']: {'claimed': False} for flow in flows}
        )

    def test_claim_queue_items(self):
        queue_kwargs = {'label': 'initial_label',
                        'queue_spec': {'item_type': 'Flow'}}
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        flows = [
            self.dao.create_item(
                item_type='Flow',
                kwargs={'label': 'flow_%s' % i}
            )
            for i in range(3)
        ]
        claimed_flows = self._claim_flows(queue=queue)
        actual_keys = sorted([flow['key'] for flow in claimed_flows])
        expected_keys = sorted([flow['key'] for flow in flows])
        self.assertEqual(actual_keys, expected_keys)

    def test_locks(self):
        queue_kwargs = {'label': 'initial_label',
                        'queue_spec': {'item_type': 'Flow'}}
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        unlocked_flow = self.dao.create_item(item_type='Flow',
                                             kwargs={'label': 'unlocked'})
        flow_to_unlock = self.dao.create_item(item_type='Flow',
                                              kwargs={'label': 'to_unlock',
                                                      'num_running_tasks': 1})
        locker_key = 'the_locker'
        self.dao.create_item(item_type='Lock',
                             kwargs={'lockee_key': flow_to_unlock['key'],
                                     'locker_key': locker_key})
        claimed_flows = self._claim_flows(queue=queue)
        self.assert_flow_lists_match(claimed_flows, [unlocked_flow])
        self._release_flows(flows=claimed_flows)

        self.dao.release_locks(locker_keys=[locker_key])
        claimed_flows = self._claim_flows(queue=queue)
        self.assert_flow_lists_match(claimed_flows,
                                     [unlocked_flow, flow_to_unlock])
        self._release_flows(flows=claimed_flows)

    def assert_flow_lists_match(self, flows_a, flows_b):
        self.assertEqual(sorted([flow['key'] for flow in flows_a]),
                         sorted([flow['key'] for flow in flows_b]))
