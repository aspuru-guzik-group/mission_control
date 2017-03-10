import unittest

class BaseTestCase(unittest.TestCase): pass

class ParserTestCase(BaseTestCase):
    def test_generates_expected_chemthings(self):
        input_dir = None
        expected_chemthings = [
            {
                'cml': '',
                'props': {}
            },
        ]
        expected_output_dir = None
        self.fail()

