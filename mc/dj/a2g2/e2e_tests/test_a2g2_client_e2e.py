from django.conf.urls import url, include
from django.test import TestCase, override_settings

from mc_utils import test_utils

from ..a2g2_dj.models import ChemThing
from ..a2g2_dj import urls as _a2g2_dj_urls
from ..a2g2_client.a2g2_client import A2G2_Client

BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_a2g2_dj_urls.__name__)),
]

class BaseTestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client(
            request_client=self.client,
            methods_to_patch=['get', 'post', 'patch']
        )
        self.a2g2_client = A2G2_Client(base_url='/' + BASE_PATH,
                                       request_client=self.client)

@override_settings(ROOT_URLCONF=__name__)
class CreateChemThingsTestCase(BaseTestCase):
    def test_creates_chemthings(self):
        self.assertEqual(ChemThing.objects.count(), 0)
        self.a2g2_client.create_chemthing(chemthing={})
        self.assertEqual(ChemThing.objects.count(), 1)

    def test_response_has_uuid(self):
        created_chemthing = self.a2g2_client.create_chemthing(chemthing={})
        self.assertTrue(created_chemthing['uuid'] is not None)

@override_settings(ROOT_URLCONF=__name__)
class GetCountsTestCase(BaseTestCase):
    def test_gets_chemthing_counts(self):
        for i in range(3): ChemThing.objects.create()
        self.assertEqual(self.a2g2_client.get_counts().get('ChemThing'),
                         ChemThing.objects.count())
