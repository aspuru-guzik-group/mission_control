import os
import subprocess
import textwrap
import unittest


class TestExample(unittest.TestCase):
    def test_example(self):
        cmd = 'cd {this_dir} && python entrypoint.py'.format(
            this_dir=os.path.dirname(os.path.abspath(__file__))
        )
        stdout = subprocess.check_output(cmd, shell=True).decode()
        expected_stdout = textwrap.dedent(
            '''
            I am task_1.
            creating job_record
            getting job_record, job_meta: {'key': 'some_key'}
            '''
        ).lstrip()
        self.assertEqual(stdout, expected_stdout)
