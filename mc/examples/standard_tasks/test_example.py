import os
import subprocess
import unittest


class TestExample(unittest.TestCase):
    def test_example(self):
        cmd = 'cd {this_dir} && python entrypoint.py'.format(
            this_dir=os.path.dirname(os.path.abspath(__file__))
        )
        stdout = subprocess.check_output(cmd, shell=True).decode()
        expected_stdout = 'message set by task_1\n'
        self.assertEqual(stdout, expected_stdout)

    def test_add_test_for_switch_task(self):
        self.fail()
