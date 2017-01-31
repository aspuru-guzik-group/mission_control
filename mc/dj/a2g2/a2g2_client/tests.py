from unittest import TestCase
from unittest.mock import call, patch, Mock

from . import a2g2_client

class BaseTestCase(TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.patchers = {}
        self.decorate_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}
        self.a2g2_client = a2g2_client.A2G2_Client(base_url=self.base_url)

    def decorate_patchers(self):
        self.patchers['requests'] = patch.object(a2g2_client, 'requests')

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

class CreateMolTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mol = Mock()

    def test_makes_post(self):
        self.a2g2_client.create_mol(mol=self.mol)
        expected_url = self.base_url + 'mols/'
        self.assertEqual(self.mocks['requests'].post.call_args,
                         call(expected_url, data=self.mol,
                              content_type='application/json'))

    def test_returns_response(self):
        result = self.a2g2_client.create_mol(mol=self.mol)
        self.assertEqual(
            result,
            self.mocks['requests'].post.return_value.json.return_value)

class GetCountsTestCase(BaseTestCase):
    def test_generates_expected_request(self):
        self.a2g2_client.get_counts()
        expected_url = self.base_url + 'counts/'
        self.assertEqual(self.mocks['requests'].get.call_args,
                         call(expected_url))

    def test_returns_response(self):
        result = self.a2g2_client.get_counts()
        self.assertEqual(
            result,
            self.mocks['requests'].get.return_value.json.return_value)

class QueryTestCase(BaseTestCase):
    def test_queries_mols(self):
        result = self.a2g2_client.query(q={'collection': 'mols'})
        expected_url = self.base_url + 'mols/'
        self.assertEqual(self.mocks['requests'].get.call_args,
                         call(expected_url))
        self.assertEqual(
            result,
            self.mocks['requests'].get.return_value.json.return_value)
