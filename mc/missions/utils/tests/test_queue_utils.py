import unittest
from unittest.mock import call, DEFAULT, MagicMock, patch

from .. import queue_utils as _queue_utils

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.queue = MagicMock()
        self.queue_spec = MagicMock()
        self.models = MagicMock()
        self.serializers = MagicMock()

class SerializeQueueItemsTestCase(BaseTestCase):
    def test_dispatches_to_serializer(self, *mock_args, **mock_kwargs):
        items = [MagicMock() for i in range(3)]
        result = _queue_utils.serialize_queue_items(
            queue=self.queue, items=items, serializers=self.serializers)
        expected_serializer = \
                self.serializers[self.queue.queue_spec['item_type']]
        expected_result = [expected_serializer(item).data for item in items]
        self.assertEqual(result, expected_result)

class ClaimQueueItemsTestCase(BaseTestCase):
    @patch.multiple(_queue_utils, generate_base_queryset_for_queue=DEFAULT,
                    get_claimed_filter=DEFAULT, claim_queryset_items=DEFAULT)
    def test_returns_updated_base_queryset_results(self, *mock_args,
                                                   **mock_kwargs):
        result = _queue_utils.claim_queue_items(queue=self.queue,
                                                models=self.models)
        self.assertEqual(
            _queue_utils.generate_base_queryset_for_queue.call_args,
            call(queue=self.queue, models=self.models)
        )
        expected_base_qset = \
                _queue_utils.generate_base_queryset_for_queue.return_value
        self.assertEqual(expected_base_qset.exclude.call_args,
                         call(_queue_utils.get_claimed_filter.return_value))
        expected_items_to_claim = expected_base_qset.exclude.return_value
        self.assertEqual(_queue_utils.claim_queryset_items.call_args,
                         call(queryset=expected_items_to_claim))
        self.assertEqual(result, list(expected_items_to_claim))

class GenerateBaseQuerySetForQueue(BaseTestCase):
    def test_uses_root_model(self, *mock_args, **mock_kwargs):
        result = _queue_utils.generate_base_queryset_for_queue(
            queue=self.queue, models=self.models)
        expected_root_model = self.models[self.queue.queue_spec['item_type']]
        self.assertEqual(result, expected_root_model.objects.filter())
