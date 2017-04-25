from django.test import TestCase

from unittest.mock import call, DEFAULT, MagicMock, patch

from .. import queue_utils as _queue_utils

class BaseTestCase(TestCase):
    def setUp(self):
        self.queue = MagicMock()
        self.query_params = MagicMock()
        self.queue_spec = MagicMock()

class SerializeQueueItemsTestCase(BaseTestCase):
    @patch.multiple(_queue_utils, get_serializer_for_queue=DEFAULT)
    def test_dispatches_to_serializer(self, *mock_args, **mock_kwargs):
        queue_items = [MagicMock() for i in range(3)]
        result = _queue_utils.serialize_queue_items(queue=self.queue,
                                                    queue_items=queue_items)
        expected_serializer = \
                _queue_utils.get_serializer_for_queue.return_value
        self.assertEqual(_queue_utils.get_serializer_for_queue.call_args,
                         call(queue=self.queue))
        expected_result = [expected_serializer(item).data
                           for item in queue_items]
        self.assertEqual(result, expected_result)

class GetSerializerForQueueTestCase(BaseTestCase):
    @patch.multiple(_queue_utils, get_root_model_for_queue_spec=DEFAULT,
                    missions_serializers=DEFAULT)
    def test_gets_serializer_from_registry(self, *mock_args, **mock_kwargs):
        result = _queue_utils.get_serializer_for_queue(queue=self.queue)
        expected_root_model = \
                _queue_utils.get_root_model_for_queue_spec.return_value
        self.assertEqual(_queue_utils.get_root_model_for_queue_spec.call_args,
                         call(queue_spec=self.queue.queue_spec))
        self.assertEqual(result,
                         _queue_utils.missions_serializers[expected_root_model])

class GetRootModelForQueueSpec(BaseTestCase):
    @patch.multiple(_queue_utils, _apps=DEFAULT)
    def test_dispatches_to_dj_apps(self, *mock_args, **mock_kwargs):
        result = _queue_utils.get_root_model_for_queue_spec(
            queue_spec=self.queue_spec)
        self.assertEqual(result, _queue_utils._apps.get_model.return_value)
        self.assertEqual(_queue_utils._apps.get_model.call_args,
                         call(self.queue_spec['root_model_spec']))

def ClaimQueueItemsTestCase(BaseTestCase):
    @patch.multiple(_queue_utils, generate_base_queryset_for_queue_spec=DEFAULT,
                    get_claimed_filter=DEFAULT, claim_queryset_items=DEFAULT)
    def test_returns_updated_base_queryset_results(self, *mock_args,
                                                   **mock_kwargs):
        result = _queue_utils.claim_queue_items(queue=self.queue,
                                                query_params=self.query_params)
        expected_base_qset = \
                _queue_utils.generate_queryset_for_queue_spec.return_value
        self.assertEqual(
            _queue_utils.generate_queryset_for_queue_spec.call_args,
            call(queue_spec=self.queue.queue_spec,
                 query_params=self.query_params))
        self.assertEqual(result, expected_base_qset.exclude.return_value)
        self.assertEqual(expected_base_qset.exclude.call_args,
                         call(_queue_utils.get_claimed_filter.return_value))
        self.assertEqual(_queue_utils.claim_queryset_items.call_args,
                         call(queryset=expected_base_qset))

class GenerateBaseQuerySetForQueueSpecTestCase(BaseTestCase):
    @patch.multiple(_queue_utils, get_root_model_for_queue_spec=DEFAULT)
    def test_uses_root_model(self, *mock_args, **mock_kwargs):
        result = _queue_utils.generate_base_queryset_for_queue_spec(
            queue_spec=self.queue_spec)
        expected_root_model = \
                _queue_utils.get_root_model_for_queue_spec.return_value
        self.assertEqual(result, expected_root_model.objects.filter())
        self.assertEqual(_queue_utils.get_root_model_for_queue_spec.call_args,
                         call(queue_spec=self.queue_spec))
