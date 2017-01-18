from django.conf.urls import url, include
from django.test import TestCase, override_settings

from mc_utils import test_utils

from ..a2g2_dj.models import Mol
from ..a2g2_dj import urls as _a2g2_dj_urls
from ..a2g2_client.a2g2_client import A2G2_Client

BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_a2g2_dj_urls.__name__)),
]

class BaseTestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client(request_client=self.client,
                                        json_methods=['get', 'post', 'patch'])

        self.a2g2_client = A2G2_Client(base_url='/' + BASE_PATH,
                                       request_client=self.client)

@override_settings(ROOT_URLCONF=__name__)
class CreateMolsTestCase(BaseTestCase):
    def test_creates_mols(self):
        self.assertEqual(Mol.objects.count(), 0)
        self.a2g2_client.create_mol(mol={})
        self.assertEqual(Mol.objects.count(), 1)

    def test_response_has_uuid(self):
        created_mol = self.a2g2_client.create_mol(mol={})
        self.assertTrue(created_mol['uuid'] is not None)

@override_settings(ROOT_URLCONF=__name__)
class GetCountsTestCase(BaseTestCase):
    def test_gets_mol_counts(self):
        for i in range(3): Mol.objects.create()
        self.assertEqual(self.a2g2_client.get_counts().get('Mol'),
                         Mol.objects.count())
