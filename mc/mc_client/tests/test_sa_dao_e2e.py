import json
import unittest

from mc.mc_client import sa_dao


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
        self.job_kwargs = {'label': 'initial_label'}

    def _create_job(self):
        return self.dao.create_item(item_type='Job', kwargs=self.job_kwargs)

    def test_create_job(self):
        created_job = self._create_job()
        self.assertEqual(created_job['label'], self.job_kwargs['label'])
        self.assertTrue(created_job['key'] is not None)

    def test_patch_job(self):
        created_job = self._create_job()
        patched_job = self.dao.patch_item(
            item_type='Job',
            key=created_job['key'],
            patches={'label': 'label2'}
        )
        self.assertEqual(patched_job['label'], 'label2')

class QueueTestCase(BaseTestCase):
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
        queue_kwargs = {
            'label': 'initial_label',
            'queue_spec': json.dumps({'item_type': 'Flow'})
        }
        queue = self._create_queue(queue_kwargs=queue_kwargs)
        flows = [
            self.dao.create_item(
                item_type='Flow',
                kwargs={'label': 'flow_%s' % i}
            )
            for i in range(3)
        ]
        claimed_flows = \
                self.dao.claim_queue_items(queue_key=queue['key'])['items']
        actual_keys = sorted([flow['key'] for flow in claimed_flows])
        expected_keys = sorted([flow['key'] for flow in flows])
        self.assertEqual(actual_keys, expected_keys)
