import unittest

from ..job_spec_client import JobSpecClient


class JobSpecClientTestCase(unittest.TestCase):
    def setUp(self):
        self.job_spec_client = JobSpecClient()

    def test_has_required_methods(self):
        required_methods = [
            'claim_job_specs',
            'fetch_job_specs',
            'update_job_specs'
        ]
        method_existences = {
            method_name: getattr(
                self.job_spec_client, method_name, None) is not None
            for method_name in required_methods
        }
        expected_method_existences = {
            method_name: True
            for method_name in required_methods
        }
        self.assertEqual(method_existences, expected_method_existences)

if __name__ == '__main__':
    unittest.main()
