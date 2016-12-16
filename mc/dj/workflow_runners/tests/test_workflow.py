import unittest
#from unittest.mock import call, DEFAULT, MagicMock, patch

from ..workflow import Workflow


class BaseTestCase(unittest.TestCase):
    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start()
                 for key, patcher in patchers.items()}
        return mocks

class AddNodesTestCase(BaseTestCase):
    def test_adds_nodes(self):
        self.fail()

class AddNodeTestCase(BaseTestCase):
    def test_creates_default_id(self):
        self.fail()

class AddEdgesTestCase(BaseTestCase):
    def test_adds_edges(self):
        self.fail()

class AddEdgeTestCase(BaseTestCase):
    def test_add_edge(self):
        self.fail()

if __name__ == '__main__':
    unittest.main()
