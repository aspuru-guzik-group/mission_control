import unittest

from ..nodes.b3lyp import B3LYPNode
from ..nodes.job_wrapper import JobWrapperNode


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()

    def generate_node(self, **node_kwargs):
        return B3LYPNode(**node_kwargs)

class InheritanceTestCase(BaseTestCase):
    def test_inherits_from_JobWrapper(self):
        self.assertTrue(issubclass(B3LYPNode, JobWrapperNode))

class GetJobInputTestCase(BaseTestCase):
    def test_generates_expected_job_input(self):
        node = self.generate_node(data={
            'input': {
                'xyz': 'some xyz',
                'b3lyp_spec': 'some b3lyp_spec',
            }
        })
        expected_job_input = {
            'job_type': 'b3lyp',
            'job_spec': {
                'xyz': node.data['input']['xyz'],
                'b3lyp': node.data['input']['b3lyp_spec']
            }
        }
        self.assertEqual(node.get_job_input(), expected_job_input)

if __name__ == '__main__':
    unittest.main()
