import json

from django.conf.urls import url, include
from django.test import TestCase
from django.test.utils import override_settings
from rest_framework.test import APITestCase

from ..models import Mol
from ..serializers import MolSerializer
from .. import urls as _urls

BASE_PATH = 'test_api/'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_urls.__name__))
]

@override_settings(ROOT_URLCONF=__name__)
class ListMolsTestCase(APITestCase):
    def test_list_mols(self):
        mols = [Mol.objects.create() for i in range(3)]
        response = self.client.get('/' + BASE_PATH + 'mols/')
        expected_data = [MolSerializer(mol).data for mol in mols]
        self.assertEqual(sorted(response.data, key=lambda j: j['uuid']),
                         sorted(expected_data, key=lambda j:j['uuid']))

@override_settings(ROOT_URLCONF=__name__)
class PatchMolTestCase(APITestCase):
    def setUp(self):
        self.mols = [Mol.objects.create() for i in range(1)]

    def test_patch_mol(self):
        mol_to_patch = self.mols[0]
        new_values = {'props': {'k1': 'v1'}}
        response = self.client.patch(
            '/' + BASE_PATH + 'mols/%s/' % mol_to_patch.uuid, new_values,
            format='json')
        self.assertEqual(response.status_code, 200)
        patched_mol = Mol.objects.get(uuid=mol_to_patch.uuid)
        patched_mol_attrs = {attr: getattr(patched_mol, attr)
                             for attr in new_values.keys()}
        self.assertEqual(patched_mol_attrs, new_values)

@override_settings(ROOT_URLCONF=__name__)
class GetCountsTestCase(TestCase):
    def setUp(self):
        self.create_mols()

    def create_mols(self):
        for i in range(3): Mol.objects.create()

    def _get_counts(self):
        return self.client.get('/' + BASE_PATH + 'counts/',
                               content_type='application/json')

    def test_response_data(self):
        response = self._get_counts()
        self.assertEqual(response.status_code, 200)
        expected_data = {'Mol': Mol.objects.count()}
        self.assertEqual(json.loads(response.content.decode()),
                         expected_data)
