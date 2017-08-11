from collections import defaultdict
import unittest
from unittest.mock import call, MagicMock

from .. import base_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.item_type = MagicMock()

        class MyDao(base_dao.BaseDao):
            def delete_items(self, *args, **kwargs): pass
            def get_flow_queue_items_to_claim(self, *args, **kwargs): pass
            def get_items(self, *args, **kwargs): pass
            def patch_item(self, *args, **kwargs): pass

        self.dao = MyDao()

    def setup_dao_mocks(self, attrs=None):
        for attr in attrs: setattr(self.dao, attr, MagicMock())

    def generate_mocks(self, n=3): return [MagicMock() for i in range(n)]

class ClaimQueueItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_key = MagicMock()
        self.setup_dao_mocks(attrs=['get_item_by_key',
                                    'get_queue_items_to_claim', 'patch_items'])
        self.expected_queue  = self.dao.get_item_by_key.return_value
        self.dao.get_queue_items_to_claim.return_value = self.generate_mocks()
        self.result = self.dao.claim_queue_items(queue_key=self.queue_key)

    def test_gets_queue(self):
        self.assertEqual(self.dao.get_item_by_key.call_args,
                         call(item_type='queue', key=self.queue_key))

    def test_gets_items_to_claim(self):
        self.assertEqual(self.dao.get_queue_items_to_claim.call_args,
                         call(queue=self.dao.get_item_by_key.return_value))

    def test_patches_claimed_items(self):
        expected_items = self.dao.get_queue_items_to_claim.return_value
        self.assertEqual(
            self.dao.patch_items.call_args,
            call(item_type=self.expected_queue['queue_spec']['item_type'],
                 keyed_patches={item['key']: {'claimed': True}
                                for item in expected_items}
                )
        )

    def test_returns_claimed_items(self):
        self.assertEqual(self.result,
                         {'items': self.dao.patch_items.return_value.values()})

class GetQueueItemsToClaimTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_spec = defaultdict(MagicMock)
        self.queue = defaultdict(MagicMock, {'queue_spec': self.queue_spec})
        self.setup_dao_mocks(attrs=['get_flow_queue_items_to_claim',
                                    'default_get_queue_items_to_claim'])

    def _get_queue_items_to_claim(self):
        return self.dao.get_queue_items_to_claim(queue=self.queue)

    def test_dispatches_to_flow_handler_for_flow_queues(self):
        self.queue_spec['item_type'] = 'flow'
        result = self._get_queue_items_to_claim()
        self.assertEqual(self.dao.get_flow_queue_items_to_claim.call_args,
                         call(queue=self.queue))
        self.assertEqual(result,
                         self.dao.get_flow_queue_items_to_claim.return_value)

    def test_dispatches_to_default_handler_for_non_flow_queues(self):
        result = self._get_queue_items_to_claim()
        self.assertEqual(self.dao.default_get_queue_items_to_claim.call_args,
                         call(queue=self.queue))
        self.assertEqual(result,
                         self.dao.default_get_queue_items_to_claim.return_value)

class DefaultGetQueueItemsToClaim(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_dao_mocks(attrs=['get_items'])
        self.queue = MagicMock()
        self.result = \
                self.dao.default_get_queue_items_to_claim(queue=self.queue)

    def test_gets_items_to_claim(self):
        self.assertEqual(
            self.dao.get_items.call_args,
            call(item_type=self.queue['queue_spec']['item_type'],
                 query={
                    'filters': [
                        {'prop': 'claimed', 'op': '=', 'arg': False},
                        {'prop': 'status', 'op': 'IN',
                         'arg': ['PENDING', 'RUNNING']}
                    ]
                 }
                )
        )
        self.assertEqual(self.result, self.dao.get_items.return_value)

if __name__ == '__main__':
    unittest.main()


