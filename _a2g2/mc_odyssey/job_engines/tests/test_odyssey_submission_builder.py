import unittest


class BaseTestCase(unittest.TestCase): pass

class BuildSubmissionTestCase(BaseTestCase):
    def test_something(self):
        self.fail("implement!")
