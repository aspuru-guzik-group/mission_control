import unittest

from ..nodes.confgen import Confgen_Node
from ..nodes.job_wrapper import JobWrapperNode


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()

    def generate_node(self, **node_kwargs):
        return Confgen_Node(**node_kwargs)

class InheritanceTestCase(BaseTestCase):
    def test_inherits_from_JobWrapper(self):
        self.assertTrue(issubclass(Confgen_Node, JobWrapperNode))

class GetJobInputTestCase(BaseTestCase):
    def test_generates_expected_job_input(self):
        node = self.generate_node(data={
            'input': {
                'smiles': 'some smiles',
                'confgen_spec': 'some confgen_spec',
            }
        })
        expected_job_input = {
            'job_type': 'confgen',
            'job_spec': {
                'smiles': node.data['input']['smiles'],
                'confgen': node.data['input']['confgen_spec']
            }
        }
        self.assertEqual(node.get_job_input(), expected_job_input)


if __name__ == '__main__':
    unittest.main()
