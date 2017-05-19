import unittest

from .. import dj_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.dao = dj_dao.DjDao(db_cfg={
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': ':memory:?foo'
        })

class FlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_kwargs = {'label': 'initial_label'}

    def _create_flow(self):
        return self.dao.create_flow(flow_kwargs=self.flow_kwargs)

    def test_create_flow(self):
        created_flow = self._create_flow()
        self.assertEqual(created_flow['label'], self.flow_kwargs['label'])
        self.assertTrue(created_flow['uuid'] is not None)

    def test_patch_flow(self):
        created_flow = self._create_flow()
        patched_flow = self.dao.patch_flow(key=created_flow['uuid'],
                                           patches={'label': 'label2'})
        self.assertEqual(patched_flow['label'], 'label2')

class JobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_kwargs = {'label': 'initial_label'}

    def _create_job(self):
        return self.dao.create_job(job_kwargs=self.job_kwargs)

    def test_create_job(self):
        created_job = self._create_job()
        self.assertEqual(created_job['label'], self.job_kwargs['label'])
        self.assertTrue(created_job['uuid'] is not None)

    def test_patch_job(self):
        created_job = self._create_job()
        patched_job = self.dao.patch_job(key=created_job['uuid'],
                                           patches={'label': 'label2'})
        self.assertEqual(patched_job['label'], 'label2')

class QueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_kwargs = {'label': 'initial_label'}

    def _create_queue(self, queue_kwargs=None):
        queue_kwargs = queue_kwargs or self.queue_kwargs
        return self.dao.create_queue(queue_kwargs=queue_kwargs)

    def test_create_queue(self):
        created_queue = self._create_queue()
        self.assertEqual(created_queue['label'], self.queue_kwargs['label'])
        self.assertTrue(created_queue['uuid'] is not None)

    def test_claim_queue_items(self):
        queue_kwargs = {
            'label': 'initial_label',
            'queue_spec': {
                'item_type': 'Flow'
            }
        }
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        flows = [self.dao.create_flow(flow_kwargs={'label': 'flow_%s' % i})
                 for i in range(3)]
        claimed_flows = self.dao.claim_queue_items(queue_key=queue['uuid'])
        actual_uuids = sorted([flow['uuid'] for flow in claimed_flows['items']])
        expected_uuids = sorted([flow['uuid'] for flow in flows])
        self.assertEqual(actual_uuids, expected_uuids)
