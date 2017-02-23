import json

from django.conf.urls import url, include
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..models import a2g2_dj_models, ChemThing
from ..serializers import ChemThingSerializer
from .. import urls as _urls

BASE_PATH = 'test_api/'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_urls.__name__))
]

@override_settings(ROOT_URLCONF=__name__)
class ListChemThingsTestCase(APITestCase):
    def test_list_chemthings(self):
        chemthings = [ChemThing.objects.create() for i in range(3)]
        response = self.client.get('/' + BASE_PATH + 'chemthings/')
        expected_data = [ChemThingSerializer(chemthing).data for chemthing in chemthings]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

@override_settings(ROOT_URLCONF=__name__)
class PatchChemThingTestCase(APITestCase):
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

@override_settings(ROOT_URLCONF=__name__)
class GetCountsTestCase(TestCase):
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

@override_settings(ROOT_URLCONF=__name__)
class FlushTestCase(TestCase):
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
