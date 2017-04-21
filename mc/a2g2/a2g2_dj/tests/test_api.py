import json
from unittest.mock import patch

from django.conf.urls import url, include
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from .. import views
from ..models import a2g2_dj_models, ChemThing
from ..serializers import ChemThingSerializer
from .. import urls as _urls

BASE_PATH = 'test_api/'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_urls.__name__))
]

@override_settings(ROOT_URLCONF=__name__)
class BaseAPITestCase(APITestCase): pass

class ListChemThingsTestCase(BaseAPITestCase):
    def test_list_chemthings(self):
        chemthings = [ChemThing.objects.create() for i in range(3)]
        response = self.client.get('/' + BASE_PATH + 'chemthings/')
        expected_data = [ChemThingSerializer(chemthing).data
                         for chemthing in chemthings]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

class QueryChemThingsTestCase(BaseAPITestCase):
    def _get_sorted_chemthings(self, chemthings=None):
        return sorted(chemthings, key=lambda c: c['uuid'])

    def run_query_test(self, wanted_chemthings=None, query_params=None):
        response = self.client.get('/' + BASE_PATH + 'chemthings/')
        expected_data = [ChemThingSerializer(chemthing).data
                         for chemthing in wanted_chemthings]
        self.assertEqual(self._get_sorted_chemthings(response.data),
                         self._get_sorted_chemthings(expected_data))

    def test_filter_by_uuid(self):
        uuid = 'some_uuid'
        self.run_query_test(
            wanted_chemthings=[ChemThing.objects.create(uuid=uuid)],
            query_params={'uuid': uuid}
        )

    def test_filter_by_tag(self):
        tags = ['tag_%s' % i for i in range(2)]
        tagged_chemthings = self.generate_tagged_chemthings(
            tags=tags, tags_per_chemthing=2)
        chemthings_by_tag = {
            tag: [chemthing for chemthing in tagged_chemthings
                  if tag in chemthing.tags]
            for tag in tags
        }
        for tag in tags:
            wanted_chemthings = chemthings_by_tag[tag]
            self.run_query_test(wanted_chemthings=wanted_chemthings,
                                query_params={'tag': tag})

    def generate_tagged_chemthings(self, tags=None, tags_per_chemthing=1,
                                   num_chemthings=None):
        num_chemthings = num_chemthings or (2 * len(tags) + 1)
        chemthings = []
        for i in range(num_chemthings):
            tags = [tags[(i + j) % len(tags)]
                    for j in range(tags_per_chemthing)]
            chemthings.append(ChemThing.objects.create(tags=tags))
        return chemthings

class PatchChemThingTestCase(BaseAPITestCase):
    def setUp(self):
        self.chemthings = [ChemThing.objects.create() for i in range(1)]

    def test_patch_chemthing(self):
        chemthing_to_patch = self.chemthings[0]
        new_values = {'props': {'k1': 'v1'}}
        response = self.client.patch(
            '/' + BASE_PATH + 'chemthings/%s/' % chemthing_to_patch.uuid, new_values,
            format='json')
        self.assertEqual(response.status_code, 200)
        patched_chemthing = ChemThing.objects.get(uuid=chemthing_to_patch.uuid)
        patched_chemthing_attrs = {attr: getattr(patched_chemthing, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_chemthing_attrs, new_values)

class PostChemThingBulkActionsTestCase(BaseAPITestCase):
    @patch.object(views, '_a2g2_dj_utils')
    def test_list_chemthings(self, mock_a2g2_dj_utils):
        mock_result = 'mock_result'
        mock_a2g2_dj_utils.process_serialized_chemthing_actions\
                .return_value = mock_result
        serialized_bulk_actions = 'mock_serialized_bulk_actions'
        bulk_url = '/' + BASE_PATH + 'chemthings/_bulk/'
        response = self.client.generic('POST', bulk_url,
                                       serialized_bulk_actions)
        self.assertEqual(response.json(), mock_result)

class GetCountsTestCase(BaseAPITestCase):
    def setUp(self):
        self.create_chemthings()

    def create_chemthings(self):
        for i in range(3): ChemThing.objects.create()

    def _get_counts(self):
        return self.client.get('/' + BASE_PATH + 'counts/')

    def test_response_data(self):
        response = self._get_counts()
        self.assertEqual(response.status_code, 200)
        expected_data = {'ChemThing': ChemThing.objects.count()}
        self.assertEqual(json.loads(response.content.decode()),
                         expected_data)

class FlushTestCase(BaseAPITestCase):
    def setUp(self):
        self.create_models()

    def create_models(self):
        for model in a2g2_dj_models:
            for i in range(3): model.objects.create()

    def _flush(self):
        return self.client.get('/' + BASE_PATH + 'flush/')

    def test_response_data(self):
        response = self._flush()
        self.assertEqual(response.status_code, 200)
        for model in a2g2_dj_models:
            self.assertEqual(model.objects.count(), 0)
        expected_flush_results = {model.__name__: 'flushed'
                                  for model in a2g2_dj_models}
        self.assertEqual(json.loads(response.content.decode()), 
                         expected_flush_results)
