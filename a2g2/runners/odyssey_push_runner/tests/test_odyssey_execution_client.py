import unittest

class BaseTestCase(unittest.TestCase):
    pass

@unittest.skip("fix later")
class DecorateExecutionStateTestCase(BaseTestCase):
    def test_gets_stdout(self):
        self.fail()
