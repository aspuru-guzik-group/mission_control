import unittest

class BaseTestCase(unittest.TestCase):
    pass

@unittest.skip("fix later")
class GenerateTaskHandler(BaseTestCase):
    def test_dispatches_to_a2g2_task_handler_compile(self):
        self.fail()

@unittest.skip("fix later")
class ExposeOuputsTestCase(BaseTestCase):
    def test_exposes_artifact(self):
        self.fail()

    def test_exposes_stdout(self):
        self.fail()
