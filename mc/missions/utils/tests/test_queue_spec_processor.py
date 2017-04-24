from django.test import TestCase

from unittest.mock import call, MagicMock, patch

from .. import queue_spec_processor


class BaseTestCase(TestCase):
    def setUp(self):
        self.processor = queue_spec_processor.QueueSpecProcessor()
        self.query_spec = MagicMock()

class GenerateQuerySetTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.processor.generate_base_queryset = MagicMock()
        self.expected_base_queryset = \
                self.processor.generate_base_queryset().return_value

    def test_generates_expected_queryset(self):
        queryset = self.processor.generate_queryset()
        self.assertEqual(queryset, self.expected_base_queryset)

class GenerateBaseQuerySetTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.processor.get_model = MagicMock()

    def _generate_base_queryset(self):
        self.base_queryset = self.processor.generate_base_queryset(
            query_spec=self.query_spec)

    def test_uses_root_model(self):
        self._generate_base_queryset()
        self.assertEqual(self.base_queryset,
                         self.processor.get_model.return_value.filter)
        self.assertEqual(self.processor.get_model.call_args,
                         call(model_spec=self.query_spec['root_model_spec']))

class GetModel(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.model_spec = MagicMock()

    @patch.object(queue_spec_processor, '_apps')
    def test_dispatches_to_apps(self, _apps):
        self.assertEqual(self.processor.get_model(model_spec=self.model_spec),
                         _apps.get_model.return_value)
        self.assertEqual(_apps.get_model.call_args, self.model_spec)
