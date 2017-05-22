import unittest
from unittest.mock import call, MagicMock

from .. import dj_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.item_type = MagicMock()
        self.dao = dj_dao.DjDao(
            models=MagicMock(),
            serializers=MagicMock()
        )

    def generate_mocks(self, n=3): return [MagicMock() for i in range(n)]

class CreateItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_item_model_cls = MagicMock()
        self.dao.serialize_item_models = MagicMock()
        self.kwargs = MagicMock()
        self.result = self.dao.create_item(item_type=self.item_type,
                                           kwargs=self.kwargs)

    def test_gets_model_cls(self):
        self.assertEqual(self.dao.get_item_model_cls.call_args,
                         call(item_type=self.item_type))

    def test_creates_model_using_db_id(self):
        expected_model_cls = self.dao.get_item_model_cls.return_value
        self.assertEqual(expected_model_cls.objects.db_manager.call_args,
                         call(self.dao.db_id))
        self.assertEqual(
            expected_model_cls.objects.db_manager.return_value.create.call_args,
            call(**self.kwargs)
        )

    def test_returns_serialized_model(self):
        expected_model = self.dao.get_item_model_cls.return_value\
                .objects.db_manager.return_value.create.return_value
        self.assertEqual(
            self.dao.serialize_item_models.call_args,
            call(item_type=self.item_type, item_models=[expected_model])
        )
        self.assertEqual(self.result,
                         self.dao.serialize_item_models.return_value[0])

class GetItemModelClsTestCase(BaseTestCase):
    def test_gets_model_cls_from_dict(self):
        result = self.dao.get_item_model_cls(item_type=self.item_type)
        self.assertEqual(result, self.dao.models[self.item_type])

class SerializeItemModels(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_item_serializer = MagicMock()
        self.item_models = self.generate_mocks()
        self.result = self.dao.serialize_item_models(
            item_type=self.item_type, item_models=self.item_models)

    def test_gets_serializer(self):
        self.assertEqual(self.dao.get_item_serializer.call_args,
                         call(item_type=self.item_type))

    def test_returns_serialized_items(self):
        expected_serializer = self.dao.get_item_serializer.return_value
        expected_result = [expected_serializer(item_model).data
                           for item_model in self.item_models]
        self.assertEqual(self.result, expected_result)

class GetItemSerializer(BaseTestCase):
    def test_gets_serializer_from_dict(self):
        result = self.dao.get_item_serializer(item_type=self.item_type)
        self.assertEqual(result, self.dao.serializers[self.item_type])

class GetItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.get_item_models = MagicMock()
        self.dao.serialize_item_models = MagicMock()
        self.result = self.dao.get_items(item_type=self.item_type,
                                         query=self.query)

    def test_gets_item_models(self):
        self.assertEqual(self.dao.get_item_models.call_args,
                         call(item_type=self.item_type, query=self.query))

    def test_returns_serialized_item_models(self):
        self.assertEqual(
            self.dao.serialize_item_models.call_args,
            call(item_type=self.item_type,
                 item_models=self.dao.get_item_models.return_value)
        )
        self.assertEqual(self.result,
                         self.dao.serialize_item_models.return_value)

class GetItemModelsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.item_query_to_queryset = MagicMock()
        self.result = self.dao.get_item_models(item_type=self.item_type,
                                               query=self.query)

    def test_gets_qset(self):
        self.assertEqual(self.dao.item_query_to_queryset.call_args,
                         call(item_type=self.item_type, query=self.query))

    def test_returns_qset_result(self):
        expected_qset = self.dao.item_query_to_queryset.return_value
        self.assertEqual(self.result, expected_qset.all())

class ItemQueryToQuerysetTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.validate_query = MagicMock()
        self.dao.get_item_model_cls = MagicMock()
        self.dao.get_dj_filter_kwargs_for_query = MagicMock()
        self.result = self.dao.item_query_to_queryset(item_type=self.item_type,
                                                      query=self.query)

    def test_validates_query(self):
        self.assertEqual(self.dao.validate_query.call_args,
                         call(query=self.query))

    def test_returns_queryset(self):
        self.assertEqual(self.dao.get_item_model_cls.call_args,
                         call(item_type=self.item_type))
        expected_item_model_cls = self.dao.get_item_model_cls.return_value
        self.assertEqual(expected_item_model_cls.objects.using.call_args,
                         call(self.dao.db_id))
        expected_base_qs = expected_item_model_cls.objects.using.return_value
        self.assertEqual(self.dao.get_dj_filter_kwargs_for_query.call_args,
                         call(query=self.query))
        expected_filter_kwargs = \
                self.dao.get_dj_filter_kwargs_for_query.return_value
        self.assertEqual(expected_base_qs.filter.call_args,
                         call(**expected_filter_kwargs))
        self.assertEqual(self.result, expected_base_qs.filter.return_value)

class PatchItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.key = MagicMock()
        self.patches = MagicMock()
        self.dao.get_item_model_by_key = MagicMock()
        self.dao.patch_model = MagicMock()
        self.dao.serialize_item_models = MagicMock()
        self.result = self.dao.patch_item(item_type=self.item_type,
                                          key=self.key, patches=self.patches)


    def test_gets_model(self):
        self.assertEqual(self.dao.get_item_model_by_key.call_args,
                         call(item_type=self.item_type, key=self.key))

    def test_patches_model(self):
        self.assertEqual(
            self.dao.patch_model.call_args,
            call(model=self.dao.get_item_model_by_key.return_value,
                 patches=self.patches)
        )

    def test_returns_serialized_patched_model(self):
        self.assertEqual(
            self.dao.serialize_item_models.call_args,
            call(item_type=self.item_type,
                 item_models=[self.dao.patch_model.return_value])
        )
        self.assertEqual(self.result,
                         self.dao.serialize_item_models.return_value[0])

class GetItemModelByKeyTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_item_models = MagicMock()
        self.key = MagicMock()
        self.result = self.dao.get_item_model_by_key(item_type=self.item_type,
                                                     key=self.key)

    def test_dispatches_to_get_item_models(self):
        self.assertEqual(
            self.dao.get_item_models.call_args,
            call(item_type=self.item_type,
                 query={
                     'filters': [
                         {'field': 'key', 'operator': '=', 'value': self.key}
                     ]
                 }
                )
        )
        self.assertEqual(self.result, self.dao.get_item_models.return_value[0])

class FlushTestCase(BaseTestCase):
    def test_dispatches_to_model_deletes(self):
        self.dao.flush_mc_db()
        for item_type in self.dao.ITEM_TYPES:
            expected_model_cls = \
                    self.dao.get_item_model_cls(item_type=item_type)
            self.assertEqual(
                expected_model_cls.objects.using.call_args,
                call(self.dao.db_id)
            )
            self.assertEqual(
                expected_model_cls.objects.using.return_value.all().delete\
                .call_args,
                call()
            )

class ClaimQueueItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_key = MagicMock()
        self.dao.get_item_model_by_key = MagicMock()
        self.dao.get_items = MagicMock(return_value=self.generate_mocks())
        self.dao.patch_items = MagicMock()
        self.expected_queue = self.dao.get_item_model_by_key.return_value
        self.expected_items = self.dao.get_items.return_value
        self.result = self.dao.claim_queue_items(queue_key=self.queue_key)

    def test_gets_queue_model(self):
        self.assertEqual(self.dao.get_item_model_by_key.call_args,
                         call(item_type='Queue', key=self.queue_key))

    def test_gets_items_to_claim(self):
        self.assertEqual(
            self.dao.get_items.call_args,
            call(item_type=self.expected_queue.queue_spec['item_type'],
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
            call(item_type=self.expected_queue.queue_spec['item_type'],
                 keyed_patches={item['key']: {'claimed': True}
                                for item in self.expected_items}
                )
        )

    def test_returns_claimed_items(self):
        self.assertEqual(self.result, self.dao.patch_items.return_value)

if __name__ == '__main__':
    unittest.main()


