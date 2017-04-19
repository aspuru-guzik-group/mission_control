from unittest import TestCase
from unittest.mock import call, patch, MagicMock

from . import a2g2_client

class BaseTestCase(TestCase):
    def setUp(self):
        self.base_url = '/base/url/'
        self.request_client = MagicMock()
        self.a2g2_client = a2g2_client.A2G2_Client(
            base_url=self.base_url,
            request_client=self.request_client
        )

class CreateChemThingTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.chemthing = MagicMock()

    def test_makes_post(self):
        self.a2g2_client.create_chemthing(chemthing=self.chemthing)
        expected_url = self.base_url + 'chemthings/'
        self.assertEqual(self.request_client.post.call_args,
                         call(expected_url, json=self.chemthing))

    def test_returns_response(self):
        result = self.a2g2_client.create_chemthing(chemthing=self.chemthing)
        self.assertEqual(result,
                         self.request_client.post.return_value.json.return_value)

class PostChemThingBulkActionsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.bulk_actions = MagicMock()
        self.a2g2_client.serialize_bulk_actions = MagicMock()
        self.result = self.a2g2_client.post_chemthing_bulk_actions(
            bulk_actions=self.bulk_actions)

    def test_makes_post(self):
        expected_url = self.base_url + 'chemthings/' + '_bulk/'
        expected_data = self.a2g2_client.serialize_bulk_actions.return_value
        self.assertEqual(self.request_client.post.call_args,
                         call(expected_url, data=expected_data))

    def test_returns_response(self):
        self.assertEqual(self.result,
                         self.request_client.post.return_value.json.return_value)

class SerializeBulkActionsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.bulk_actions = list(range(3))

    @patch.object(a2g2_client, 'json')
    def test_serializes_bulk_actions(self, mock_json):
        mock_json.dumps.return_value = 'mock_json_dumps'
        serialized_bulk_actions = self.a2g2_client.serialize_bulk_actions(
            bulk_actions=self.bulk_actions)
        expected_serialization = "\n".join([mock_json.dumps(action)
                                            for action in self.bulk_actions])
        self.assertEqual(serialized_bulk_actions, expected_serialization)

class GetCountsTestCase(BaseTestCase):
    def test_generates_expected_request(self):
        self.a2g2_client.get_counts()
        expected_url = self.base_url + 'counts/'
        self.assertEqual(self.request_client.get.call_args,
                         call(expected_url))

    def test_returns_response(self):
        result = self.a2g2_client.get_counts()
        self.assertEqual(result,
                         self.request_client.get.return_value.json.return_value)

class QueryTestCase(BaseTestCase):
    def test_queries_chemthings(self):
        result = self.a2g2_client.query(q={'collection': 'chemthings'})
        expected_url = self.base_url + 'chemthings/'
        self.assertEqual(self.request_client.get.call_args, call(expected_url))
        self.assertEqual(result,
                         self.request_client.get.return_value.json.return_value)

class FlushTestCase(BaseTestCase):
    def test_gets_flush_endpoint(self):
        result = self.a2g2_client.flush_a2g2_db()
        expected_url = self.base_url + 'flush/'
        self.assertEqual(self.request_client.get.call_args, call(expected_url))
        self.assertEqual(result,
                         self.request_client.get.return_value.json.return_value)
