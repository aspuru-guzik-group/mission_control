import unittest
from unittest.mock import call, MagicMock

from .. import base_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.item_type = MagicMock()
        self.dao = base_dao.BaseDao()

    def generate_mocks(self, n=3): return [MagicMock() for i in range(n)]

class ClaimQueueItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_key = MagicMock()
        self.dao.get_queue_spec_for_queue_key = MagicMock()
        self.expected_queue_spec = \
                self.dao.get_queue_spec_for_queue_key.return_value
        self.dao.get_items = MagicMock(return_value=self.generate_mocks())
        self.expected_items = self.dao.get_items.return_value
        self.dao.patch_items = MagicMock()
        self.result = self.dao.claim_queue_items(queue_key=self.queue_key)

    def test_gets_queue_spec(self):
        self.assertEqual(self.dao.get_queue_spec_for_queue_key.call_args,
                         call(queue_key=self.queue_key))

    def test_gets_items_to_claim(self):
        self.assertEqual(
            self.dao.get_items.call_args,
            call(item_type=self.expected_queue_spec['item_type'],
                 query={
                    'filters': [
                        {'field': 'claimed', 'operator': '=', 'value': False}
                    ]
                 }
                )
        )

    def test_patches_claimed_items(self):
        self.assertEqual(
            self.dao.patch_items.call_args,
            call(item_type=self.expected_queue_spec['item_type'],
                 keyed_patches={item['key']: {'claimed': True}
                                for item in self.expected_items}
                )
        )

    def test_returns_claimed_items(self):
        self.assertEqual(self.result, self.dao.patch_items.return_value)

if __name__ == '__main__':
    unittest.main()


