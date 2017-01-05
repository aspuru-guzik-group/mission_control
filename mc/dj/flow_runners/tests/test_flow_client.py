import unittest

from ..flow_client import FlowClient


class FlowClientTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_client = FlowClient()

    def test_has_required_methods(self):
        required_methods = [
            'claim_flows',
            'fetch_claimable_flows',
            'update_flows'
        ]
        method_existences = {
            method_name: getattr(
                self.flow_client, method_name, None) is not None
            for method_name in required_methods
        }
        expected_method_existences = {
            method_name: True
            for method_name in required_methods
        }
        self.assertEqual(method_existences, expected_method_existences)

if __name__ == '__main__':
    unittest.main()
