import unittest

from ..workflow_client import WorkflowClient


class WorkflowClientTestCase(unittest.TestCase):
    def setUp(self):
        self.workflow_client = WorkflowClient()

    def test_has_required_methods(self):
        required_methods = [
            'claim_workflows',
            'fetch_claimable_workflows',
            'update_workflows'
        ]
        method_existences = {
            method_name: getattr(
                self.workflow_client, method_name, None) is not None
            for method_name in required_methods
        }
        expected_method_existences = {
            method_name: True
            for method_name in required_methods
        }
        self.assertEqual(method_existences, expected_method_existences)

if __name__ == '__main__':
    unittest.main()
