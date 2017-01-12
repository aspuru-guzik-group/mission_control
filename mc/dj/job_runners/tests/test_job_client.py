import unittest

from ..job_client import JobClient


class JobClientTestCase(unittest.TestCase):
    def setUp(self):
        self.job_client = JobClient()

    def test_has_required_methods(self):
        required_methods = [
            'claim_jobs',
            'fetch_jobs',
            'update_jobs'
        ]
        method_existences = {
            method_name: getattr(
                self.job_client, method_name, None) is not None
            for method_name in required_methods
        }
        expected_method_existences = {
            method_name: True
            for method_name in required_methods
        }
        self.assertEqual(method_existences, expected_method_existences)

if __name__ == '__main__':
    unittest.main()
