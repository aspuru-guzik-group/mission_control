import unittest

class BaseTestCase(unittest.TestCase):
    pass

class DecorateExecutionStateTestCase(BaseTestCase):
    def test_gets_stdout(self):
        self.fail()
