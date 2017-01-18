from django.test import TestCase

from mc_utils import test_utils

from .a2g2_client import A2G2_Client

class BaseTestCase(TestCase):
    def setUp(self):
        test_utils.patch_request_client(request_client=self.client)
        self.a2g2_client = A2G2_Client(request_client=self.client)

class CreateMolsTestCase(BaseTestCase):
    def test_creates_mols(self): pass

class GetCountsTestCase(BaseTestCase):
    def test_gets_mol_counts(self):
        mols = self.a2g2_client.create_mols(mols=[{} for i in range(3)])
        self.assertEqual(self.a2g2_client.get_counts().get('Mol'), len(mols))
