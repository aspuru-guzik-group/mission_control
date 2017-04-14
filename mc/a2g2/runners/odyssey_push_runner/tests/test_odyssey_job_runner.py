import unittest

class BaseTestCase(unittest.TestCase):
    pass

class ExposeOuputsTestCase(BaseTestCase):
    def test_exposes_artifact(self):
        self.fail()

    def test_exposes_stdout(self):
        self.fail()
